function fig = plot_roi_trials(roi, save_flag)
% plot measured responses for each experiment across trial types.

% check inputs
if length(roi) > 1; roi = roi(1); end
if nargin < 2; save_flag = 0; end

% get design parameters
nexps = size(roi.experiments, 2); cond_groups = group_trials(roi);
cols = lines(max(cellfun(@length, roi.model.cond_list)));

% setup figure
fig_name = [roi.nickname ' trial responses'];
fig_pos = [.1 .1 max(cellfun(@length, cond_groups)) * .2 nexps * .2];
fig = figTS(fig_name,  fig_pos);

% plot responses overlaying trials of the same type across experiments
for ee = 1:nexps
    ax(ee) = subplot(nexps, 1, ee); hold on;
    xcnt = 3; zlc = xcnt; me = []; cnt = 0; leg_str = {};
    for gg = 1:length(cond_groups)
        % get duration of trial time windows and plot custom zero line
        tl = [];
        for cc = 1:length(cond_groups{ee}{gg})
            tl(cc) = length(roi.trial_avgs{cond_groups{ee}{gg}(cc), 1, ee});
        end
        plot([zlc - 1 zlc + max(tl)], [0 0], 'k-');
        % plot measured responses in peristimulus time window
        for cc = length(cond_groups{ee}{gg}):-1:1
            x = xcnt:xcnt + tl(cc) - 1;
            ym = [roi.trial_avgs{cond_groups{ee}{gg}(cc), :, ee}]';
            col = cols(cond_groups{ee}{gg}(cc), :); cnt = cnt + 1;
            me(cnt) = lineTS(x, ym, 1, col, col, 'sem');
            leg_str(cnt) = roi.model.cond_list{ee}(cond_groups{ee}{gg}(cc));
        end
        % plot stimulus
        sy = -.5;
        for cc = 1:length(cond_groups{ee}{gg})
            stim = [xcnt + roi.model.pre_dur xcnt + tl(cc) - roi.model.post_dur];
            col = cols(cond_groups{ee}{gg}(cc), :);
            plot(stim, [sy sy], 'Color', col, 'LineWidth', 4);
            sy = sy - .2;
        end
        xcnt = xcnt + max(tl) + 3; zlc = xcnt;
    end
    conds = roi.model.cond_list{ee};
    legend(me(:), leg_str, 'Location', 'NorthWestOutside'); legend boxoff;
    title(roi.experiments{ee}); ylabel('fMRI (% signal)'); axis tight;
    set(gca, 'XColor', 'w', 'TickDir', 'out', 'YTick', -2:10, 'FontSize', 8);
end

% match y-axis limits across experiments
[xmin, xmax, ymin, ymax] = deal(0);
for ee = 1:nexps
    xlims = get(ax(ee), 'XLim'); ylims = get(ax(ee), 'YLim');
    xmin = min([xmin xlims(1)]); xmax = max([xmax xlims(2)]);
    ymin = min([ymin ylims(1)]); ymax = max([ymax ylims(2)]);
end
yticks = floor(ymin):ceil(ymax);
for ee = 1:nexps
    set(ax(ee), 'XLim', [xmin xmax], 'YLim', [ymin ymax], 'YTick', yticks);
end

% save to figures directory if applicable
if save_flag
    fpath = fullfile(roi.project_dir, 'figures');
    fname = [roi.nickname '_exps_' [roi.experiments{:}] '_' date '.fig'];
    saveas(fig, fullfile(fpath, fname), 'fig');
end

end
