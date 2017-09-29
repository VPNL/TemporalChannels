function fig = plot_roi_exps(roi, save_flag)
% plot measured responses for each trial type across experiments.

% check inputs
if length(roi) > 1; roi = roi(1); end
if nargin < 2; save_flag = 0; end

% get design parameters
nexps = size(roi.experiments, 2); cond_list = roi.model.cond_list;
all_conds = unique([cond_list{:}], 'stable'); nconds = length(all_conds);
cols = lines(nexps); cond_idxs = idx_trials(roi);

% setup figure
fig_name = [roi.nickname ' experiment responses'];
fig = figTS(fig_name, [.1 .3 .8 .4]);

% plot responses overlaying trials of the same type across experiments
xcnt = 3; zlc = xcnt;
for cc = 1:nconds
    % get duration of trial time window
    tl = length(roi.trial_avgs{cond_idxs(cc, 1), 1, 1});
    % plot custom zero line for trial
    plot([zlc - 1 zlc + tl], [0 0], 'k-');
    % plot measured response in peristimulus time window
    x = xcnt:xcnt + tl - 1;
    for ee = 1:nexps
        if cond_idxs(cc, ee) > 0
            ym = [roi.trial_avgs{cond_idxs(cc, ee), :, ee}]';
            me(ee) = lineTS(x, ym, 1, cols(ee, :), cols(ee, :), 'sem');
        end
    end
    % plot stimulus
    stim = [xcnt + roi.model.pre_dur xcnt + tl - roi.model.post_dur];
    plot(stim, [-.5 -.5], 'k-', 'LineWidth', 4);
    text(xcnt + roi.model.pre_dur - 1, -.8, all_conds{cc}, 'FontSize', 8);
    xcnt = xcnt + tl + 3; zlc = xcnt;
end

% format plot
title(roi.nickname); ylabel('fMRI (% signal)');
legend(me(:), roi.experiments(1, :)); legend boxoff; axis tight;
set(gca, 'XColor', 'w', 'TickDir', 'out', 'YTick', -2:10);

% save to figures directory if applicable
if save_flag
    fpath = fullfile(roi.project_dir, 'figures');
    fname = [roi.nickname '_exps_' [roi.experiments{:}] '_' date '.fig'];
    saveas(fig, fullfile(fpath, fname), 'fig');
end

end
