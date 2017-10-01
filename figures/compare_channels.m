function compare_channels(rois, models, fit_exps)
% Generates separate channel comparison figure for each model by plotting
% the means and SEMs of beta weights for transient vs. sustained channels.
%
% Inputs
%   1) rois: cell array of ROI names (e.g., {'IOG-faces' 'pFus-faces'})
%   2) models: cell array of model names (e.g., {'standard' '2ch-dcts'})
%   3) fit_exps: cell array of experiments used to fit the model
%
% AS 8/2017

%% check inputs
if nargin < 1 || isempty(rois)
    rois = {'IOG-faces' 'pFus-faces' 'mFus-faces' 'pSTS-faces'...
        'LOS-bodies' 'ITG-bodies' 'MTG-bodies' 'OTS-bodies' ...
        'IOS-characters' 'pOTS-characters' 'mOTS-characters'};
else
    rois = force_cell(strrep(rois, '_', '-'));
end

if nargin < 2 || isempty(models)
    models = {'2ch' '2ch-lin' '2ch-rect' '2ch-sqrt-rect' ...
        '2ch-pow' '2ch-div' '2ch-dcts' '2ch-opt'};
else
    models = force_cell(models);
end

if nargin < 3 || isempty(fit_exps)
    fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
else
    fit_exps = force_cell(fit_exps);
end

%% specify paths to results/figures directories and plotting colors
project_dir = fullfile(RAID, 'projects', 'CategoryChannels');
model_dir = fullfile(project_dir, 'model', 'TemporalChannels');
res_dir = fullfile(model_dir, 'results');
fig_dir = fullfile(model_dir, 'figures', 'channel_comparison');
cat_cols = [0 0 1; 1 0 0; .5 .5 .5];

%% load the data
betas = {};
for mm = 1:length(models)
    for rr = 1:length(rois)
        for ff = 1:size(fit_exps, 1)
            fstem = [rois{rr} '_' models{mm} '_fit' [fit_exps{ff, :}] '_*.mat'];
            d = dir(fullfile(res_dir, fstem));
            load(fullfile(d.folder, d(1).name));
            roi = roi(1); model = model(1); rps = model.run_preds;
            % norm channel predictors to have save max height
            ecs = cellfun(@isempty, rps); rps(ecs) = {zeros(1, 6)};
            mS = max(cellfun(@(X) max(max(X(:, 1:3))), rps));
            mT = max(cellfun(@(X) max(max(X(:, 4:6))), rps));
            normT = model(1).normT * mean(mS ./ mT);
            model.normT = normT;
            % refit model and get betas
            model = pred_runs(model);
            model = pred_trials(model);
            fold(ff) = tc_fit(roi, model);
        end
        roi = pool_across_folds(fold(1), fold(2));
        betas{rr, mm} = cell2mat(roi.model.betas');
    end
end

%% plot channel comparison for each stimulus category
ncats = length(roi(1).model.cat_list); pidxs = [1 5 9 4 2 6 8 10 3 7 11];
for mm = 1:length(models)
    fig = figTS([models{mm} ' channel comparison'], [.1 .1 .7 .8]);
    for rr = 1:length(rois)
        subplot(3, 4, pidxs(rr)); hold on;
        means = mean(betas{rr, mm});
        sems = std(betas{rr, mm}) / sqrt(size(betas{rr, mm}, 1) - 1);
        x_avg = means(:, ncats + 1:ncats * 2); y_avg = means(:, 1:ncats);
        x_sem = sems(:, ncats + 1:ncats * 2); y_sem = sems(:, 1:ncats);
        x_max = max(x_avg + x_sem); ymax = max(y_avg + y_sem);
        x_min = min(x_avg - x_sem); ymin = min(y_avg - y_sem);
        lims = [min([0 x_min ymin]) max([ceil(x_max) ceil(ymax)])];
        lims = [0 2];
        bl = axesTS(rois{rr}, 'Transient \beta_T', 'Sustained \beta_S', lims);
        br = range(bl); lcnt = bl(2) - br * 0.05; xl = br * 0.05;
        for cc = 1:ncats
            %oo = draw_oct(x_avg(cc), y_avg(cc), x_sem(cc), y_sem(cc), ...
            %    fill_cols(cc, :), fill_cols(cc, :)); alpha(oo, .5)
            %oo = draw_oct(xl, lcnt, br * 0.02, br * 0.02, ...
            %    fill_cols(cc, :), fill_cols(cc, :)); alpha(oo, .5);
            draw_cross(x_avg(cc), y_avg(cc), x_sem(cc), y_sem(cc), cat_cols(cc, :), 1.5);
            draw_cross(xl, lcnt, br * 0.02, br * 0.02, cat_cols(cc, :), 1);
            text(xl + br * 0.03, lcnt, roi.model.cat_list{cc}, ...
                'FontSize', 6, 'FontName', 'Helvetica', ...
                'HorizontalAlignment', 'left', ...
                'VerticalAlignment', 'middle');
            lcnt = lcnt - br * 0.06;
        end
        xlim([bl(1) 3]); ylim([bl(1) 3]);
    end
    fname = ['channel_comparison_' models{mm} '_fit' [fit_exps{:}]];
    saveas(fig, fullfile(fig_dir, [fname '.fig']), 'fig');
    print(fullfile(fig_dir, [fname '.eps']), '-depsc2', '-painters');
    close all;
end

end
