function fig = plot_roi_model(roi, save_flag)

% check inputs
if length(roi) > 1; roi = roi(1); end
if nargin < 2; save_flag = 0; end

% get design parameters, data, and predictor names
nexps = size(roi.experiments, 2); npreds = size(roi.model.betas{1}, 2);
amps = reshape([roi.model.betas{:}], npreds, [])'; TR = roi.tr;
R2 = [roi.model.varexp{:}]; R2_str = num2str(mean(R2), 3);
fit_str = [roi.model.type ' fit to ' strjoin(roi.model.fit_exps, '/')];
val_str = ['R^{2} in ' strjoin(roi.experiments, '/') ' = ' R2_str];
xlabs = label_preds(roi.model); xlabs = strrep(xlabs(1:npreds), '_', '-');
fit_str = strrep(fit_str, '_', '-'); val_str = strrep(val_str, '_', '-');

% setup figure
fig_name = [roi.nickname ' - ' roi.model.type ' model'];
fig = tch_fig(fig_name, [.1 0 .8 .2 + nexps * .15]);

% plot model weights
subplot(1 + nexps, 2, 1); hold on; tch_set_axes;
[ymin, ymax] = tch_plot_bar(amps, [0 0 0]);
axis tight; xlim([0 size(amps, 2) + 1]); ylim([ymin ymax]);
title({roi.nickname; fit_str; val_str});
xlabel('Predictor'); ylabel('Beta (% signal)');
set(gca, 'XTick', 1:npreds, 'XTickLabel', xlabs);

% plot variance explained for each session
subplot(1 + nexps, 2, 2); hold on; tch_set_axes;
[ymin, ymax] = tch_plot_bar(R2, [0 0 0]);
xlim([0 size(R2, 2) + 1]); ylim([ymin ymax]);
for ss = 1:length(roi.sessions)
    ypos = max([0 R2(ss)]) + .1; lab = num2str(R2(ss), 2);
    text(ss, ypos, lab, 'HorizontalAlignment', 'center', 'FontSize', 8);
end
title('Performance'); xlabel('Session'); ylabel('R^2'); ylim([0 1]);
set(gca, 'XTick', 1:length(roi.sessions), 'XTickLabel', roi.session_ids);

% plot measurement vs prediction for each trial type
for ee = 1:nexps
    ax(ee) = subplot(1 + nexps, 1, ee + 1); hold on; tch_set_axes;
    xcnt = 3; zlc = xcnt;
    for cc = 1:length(roi.trial_avgs(:, 1, ee))
        % plot custom zero line for trial
        tl = length(roi.trial_avgs{cc, 1, ee});
        plot([zlc - 1 zlc + tl], [0 0], 'k-');
        % plot measured response for peristimulus time window
        x = xcnt:xcnt + tl - 1; ym = [roi.trial_avgs{cc, :, ee}]';
        me = tch_plot_tc(x, ym, 1, [.9 .9 .9], [.7 .7 .7], 'std');
        % plot model prediction for peristimulus time window
        pr = tch_plot_tc(x, [roi.trial_preds_sum{cc, :, ee}]', 2, [0 0 0]);
        % plot separate channel contributions if applicable
        if roi.model.num_channels > 1
            sp = tch_plot_tc(x, [roi.trial_predsS_sum{cc, :, ee}]', 1, [0 0 1]);
            tp = tch_plot_tc(x, [roi.trial_predsT_sum{cc, :, ee}]', 1, [1 0 0]);
        end
        if roi.model.num_channels > 2
            dp = tch_plot_tc(x, [roi.trial_predsP_sum{cc, :, ee}]', 1, [0 1 0]);
        end
        % plot stimulus
        stim = [xcnt + roi.model.pre_dur / TR xcnt + tl - roi.model.post_dur / TR];
        cond_name = roi.model.cond_list{ee}(cc);
        plot(stim, [-.5 -.5], 'k-', 'LineWidth', 4);
        text(xcnt + roi.model.pre_dur / TR - 1, -1, cond_name, 'FontSize', 8);
        xcnt = xcnt + tl + 3; zlc = xcnt;
    end
    % set legend and format plot
    leg{1} = [roi.nickname ' (N = ' num2str(length(roi.sessions)) ')'];
    leg{2} = [roi.model.type ' model']; ptrs = [me pr];
    if roi.model.num_channels > 1
        leg(3:4) = {'Sustained' 'Transient'}; ptrs = [ptrs sp tp];
    end
    if roi.model.num_channels > 2
        leg{5} = 'Persistent'; ptrs = [ptrs dp];
    end
    legend(ptrs, leg, 'Location', 'NorthWestOutside'); legend boxoff;
    title([roi.experiments{1, ee}], 'FontSize', 8); ylabel('fMRI (% signal)');
    set(gca, 'XColor', 'w'); axis tight;
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
    fname = [roi.nickname '_' roi.model.type ...
        '_fit' [roi.model.fit_exps{:}] ...
        '_val' [roi.experiments{:}] ...
        '_' date '.fig'];
    saveas(fig, fullfile(fpath, fname), 'fig');
end

end
