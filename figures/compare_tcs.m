function compare_tcs(rois, exps)
% Generates model comparison figure for all ROIs by plotting trial time
% series for each condition in experiments.
% When multiple rows of fit_exps and val_exps are passed, performance is
% averaged across both splits of the data.
%
% Inputs
%   1) rois: cell array of ROI names (e.g., {'IOG-faces' 'pFus-faces'})
%   2) exps: cell array of experiments used to plot
%
% AS 8/2017

%% check inputs
if nargin < 1 || isempty(rois)
    rois = {'IOG-faces' 'pFus-faces' 'mFus-faces' 'pSTS-faces' ...
        'LOS-bodies' 'ITG-bodies' 'MTG-bodies' 'OTS-bodies' ...
        'IOS-characters' 'pOTS-characters' 'mOTS-characters'};
else
    rois = force_cell(strrep(rois, '_', '-'));
end

if nargin < 2 || isempty(exps)
    exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
    exp_names = {'Exp. 1' 'Exp. 2' 'Exp. 3'};
else
    exps = force_cell(exps);
    exp_names = exps(1, :);
end

val_exps = flipud(exps);

%% specify paths to files in results/figures directories
project_dir = fullfile(RAID, 'projects', 'CategoryChannels');
model_dir = fullfile(project_dir, 'model', 'TemporalChannels');
res_dir = fullfile(model_dir, 'results');
fig_dir = fullfile(model_dir, 'figures', 'tc_comparison');
% generate stem of data filenames
fstems = cell(1, size(exps, 1));
for ff = 1:size(exps, 1)
    fstem = ['_fit' [exps{ff, :}]];
    for ee = 1:size(val_exps, 2)
        fstem = [fstem '_val' val_exps{ff, ee}];
    end
    fstems{ff} = [fstem '.mat'];
end

%% load data for all ROIs, models, experiments, and data folds
for rr = 1:length(rois)
    for ff = 1:size(exps, 1)
        fold(ff) = load(fullfile(res_dir, [rois{rr} '_glm' fstems{ff}]));
    end
    roi(rr) = pool_across_folds(fold(1).roi(1), fold(2).roi(1));
end

%% plot tcs comparing across experiments
cond_order = [9 10 11 12 5 6 7 8 1 2 3 4];
cols = lines(length(exp_names) + 2); cols = cols(3:end, :);
scols = repmat([0 0 0; lines(2)]', 4, 1); scols = reshape(scols(:)', 3, 12)';
fig = figTS('Time series comparison across experiments', [.05 .05 .4 1.4]);
for rr = 1:length(rois)
    subplot(length(rois), 1, rr); hold on;
    % get design parameters
    nexps = size(roi(rr).experiments, 2);
    cond_list = roi(rr).model.cond_list;
    all_conds = unique([cond_list{:}], 'stable'); nconds = length(all_conds);
    all_conds = repmat({'3' '5' '10' '20'}, 1, 3);
    pre_dur = roi(rr).model.pre_dur; post_dur = roi(rr).model.post_dur;
    cond_idxs = idx_trials(roi(rr));
    % plot responses to trials of the same type across experiments
    xcnt = 3; zlc = xcnt;
    y_max = 0; y_min = -1;
    for tt = cond_order
        % get duration of trial time window
        if nexps > 1
            tl = length(roi(rr).trial_avgs{cond_idxs(tt, 1), 1, 1});
        else
            tl = length(roi(rr).trial_avgs{cond_idxs(tt, 1), 1});
        end
        sd = tl - pre_dur - post_dur;
        % plot custom zero line for trial
        plot([zlc - 1 zlc + tl], [0 0], 'k-');
        % plot measured response in peristimulus time window
        x = xcnt:xcnt + tl - 1;
        for ee = nexps:-1:1
            if cond_idxs(tt, ee) > 0
                yp_m = [roi(rr).trial_avgs{cond_idxs(tt, ee), :, ee}]';
                [me(ee), cymin, cymax] = lineTS(x, yp_m, 1, ...
                    cols(ee, :), cols(ee, :), 'sem');
                y_min = min([y_min cymin]); y_max = max([y_max cymax]);
            end
        end
        % plot stimulus
        plot([xcnt + pre_dur - 1 xcnt + tl - post_dur], [-.5 -.5], ...
            'Color', scols(tt, :), 'LineWidth', 7);
        text(xcnt + pre_dur + sd / 2, -.49, all_conds{tt}, ...
            'FontSize', 6, 'FontWeight', 'bold', 'Color', 'w', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        if sum(tt == [8 12])
            xcnt = xcnt + 6;
        end
        xcnt = xcnt + tl + 3; zlc = xcnt;
    end
    % format plot
    ylabel('fMRI (% signal)'); ylim([floor(y_min) ceil(y_max)]); ylim([-.61 3]);
    % legend(me(:), exp_names, 'Location', 'NorthEastOutside'); legend boxoff;
    title(roi(rr).nickname, 'Interpreter', 'none');
    set(gca, 'XColor', 'w', 'FontSize', 6, ...
        'TickDir', 'out', 'YTick', -1:3, 'YMinorTick', 'on');
end
fname = 'tc_comparison_across_exps';
saveas(fig, fullfile(fig_dir, [fname '.fig']), 'fig');
print(fullfile(fig_dir, [fname '.eps']), '-depsc2', '-painters');
close all;

%% plot tcs comparing across categories
cols = repmat([0 0 0; lines(2)]', 4, 1); cols = reshape(cols(:)', 3, 12)';
fig = figTS('Time series comparison across categories', [.05 .05 .4 1.4]);
for rr = 1:length(rois)
    subplot(length(rois), 1, rr); hold on;
    % get design parameters
    nexps = size(roi(rr).experiments, 2);
    cond_list = roi(rr).model.cond_list;
    all_conds = unique([cond_list{:}], 'stable'); nconds = length(all_conds);
    all_conds = repmat({'3' '5' '10' '20'}, 1, 3);
    pre_dur = roi(rr).model.pre_dur; post_dur = roi(rr).model.post_dur;
    cond_idxs = idx_trials(roi(rr));
    % plot responses to trials of the same type across experiments
    xcnt = 3; zlc = xcnt; y_max = 0; y_min = -1;
    for ee = 1:nexps
        for tt = 1:4
            % get duration of trial time window
            if nexps > 1
                tl = length(roi(rr).trial_avgs{cond_idxs(tt, 1), 1, 1});
            else
                tl = length(roi(rr).trial_avgs{cond_idxs(tt, 1), 1});
            end
            sd = tl - pre_dur - post_dur;
            % plot custom zero line for trial
            plot([zlc - 1 zlc + tl], [0 0], 'k-');
            % plot measured response in peristimulus time window
            x = xcnt:xcnt + tl - 1;
            for ii = 1:3
                ci = 4 * (ii - 1) + tt;
                if cond_idxs(ci, ee) > 0
                    yp_m = [roi(rr).trial_avgs{cond_idxs(ci, ee), :, ee}]';
                    [me(ee), cymin, cymax] = lineTS(x, yp_m, 1, ...
                        cols(ci, :), cols(ci, :), 'sem');
                    y_min = min([y_min cymin]); y_max = max([y_max cymax]);
                end
            end
            % plot stimulus
            plot([xcnt + pre_dur - 1 xcnt + tl - post_dur], [-.5 -.5], ...
                'Color', [.5 .5 .5], 'LineWidth', 7);
            text(xcnt + pre_dur + sd / 2, -.49, all_conds{tt}, ...
                'FontSize', 6, 'FontWeight', 'bold', 'Color', 'w', ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
            xcnt = xcnt + tl + 3; zlc = xcnt;
        end
        xcnt = xcnt + 6; zlc = xcnt;
    end
    % format plot
    ylabel('fMRI (% signal)'); ylim([floor(y_min) ceil(y_max)]); ylim([-.61 3]);
    % legend(me(:), exp_names, 'Location', 'NorthEastOutside'); legend boxoff;
    title(roi(rr).nickname, 'Interpreter', 'none');
    set(gca, 'XColor', 'w', 'FontSize', 6, ...
        'TickDir', 'out', 'YTick', -1:3, 'YMinorTick', 'on');
end
fname = 'tc_comparison_across_cats';
saveas(fig, fullfile(fig_dir, [fname '.fig']), 'fig');
print(fullfile(fig_dir, [fname '.eps']), '-depsc2', '-painters');
close all;

end
