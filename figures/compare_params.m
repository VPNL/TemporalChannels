function compare_params(rois, models, fit_exps, val_exps)
% Generates separate parameter comparison figure for each model by plotting
% the optimized hyperparameter values of each ROI.
% When multiple rows of fit_exps and val_exps are passed, performance is
% averaged across both splits of the data.
%
% Inputs
%   1) rois: cell array of ROI names (e.g., {'IOG-faces' 'pFus-faces'})
%   2) models: cell array of model names (e.g., {'standard' '2ch-dcts'})
%   3) fit_exps: cell array of experiments used to fit the model
%   4) val_exps: cell array of experiments to validate the model
%
% AS 9/2017

%% check inputs
if nargin < 1 || isempty(rois)
    rois = {'IOG-faces' 'pFus-faces' 'mFus-faces' 'pSTS-faces' ...
        'LOS-bodies' 'ITG-bodies' 'MTG-bodies' 'OTS-bodies' ...
        'IOS-characters' 'pOTS-characters' 'mOTS-characters'};
    roi_names = rois;
    roi_names = strrep(roi_names, '-faces', '');
    roi_names = strrep(roi_names, '-bodies', '');
    roi_names = strrep(roi_names, '-characters', '');
    rcols = [1 .5 .5; 1 .5 .5; 1 .5 .5; 1 .5 .5; ...
        .5 .5 1; .5 .5 1; .5 .5 1; .5 .5 1; ...
        .7 .7 .7; .7 .7 .7; .7 .7 .7];
    ecols = [1 0 0; 1 0 0; 1 0 0; 1 0 0; ...
        0 0 1; 0 0 1; 0 0 1; 0 0 1; ...
        0 0 0; 0 0 0; 0 0 0];
else
    rois = force_cell(rois);
    rcols = repmat([.7 .7 .7], length(rois), 1);
    ecols = repmat([0 0 0], length(rois), 1);
end

if nargin < 2 || isempty(models)
    models = {'cts-pow' 'cts-div' 'dcts' ...
        '2ch-pow' '2ch-div' '2ch-dcts' ...
        '2ch-opt' '3ch' '3ch-opt'};
else
    models = force_cell(models);
end

if nargin < 3 || isempty(fit_exps)
    fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
else
    fit_exps = force_cell(fit_exps);
end

if nargin < 4 || isempty(val_exps)
    val_exps = flipud(fit_exps);
else
    val_exps = force_cell(val_exps);
end

%% specify paths to files in results/figures directories
project_dir = fullfile(RAID, 'projects', 'CategoryChannels');
model_dir = fullfile(project_dir, 'model', 'TemporalChannels');
res_dir = fullfile(model_dir, 'results');
fig_dir = fullfile(model_dir, 'figures', 'param_comparison');
% generate stem of data filenames
for ff = 1:size(fit_exps, 1)
    fstem = ['_fit' [fit_exps{ff, :}]];
    for ee = 1:size(val_exps, 2)
        fstem = [fstem '_val' val_exps{ff, ee}];
    end
    fstems{ff} = [fstem '.mat'];
end

%% load data for all ROIs, models, experiments, and data folds
params = struct;
for mm = 1:length(models)
    for rr = 1:length(rois)
        pvs = [];
        for ff = 1:size(fit_exps, 1)
            load(fullfile(res_dir, [rois{rr} '_' models{mm} fstems{ff}]));
            if rr == 1
                params(mm).names = fieldnames(roi(1).model.params)';
            end
            for pp = 1:length(params(mm).names)
                pvs(:, pp, ff) = [roi(1).model.params.(params(mm).names{pp}){:}]';
            end
        end
        pms = mean(pvs, 3);
        for pp = 1:length(params(mm).names)
            params(mm).values{rr, pp} = pms(:, pp);
        end
    end
end

%% plot parameter comparison for each experiment and ROI
for mm = 1:length(models)
    fig = figTS([models{mm} ' parameter comparison'], [.1 .1 .3 .6]);
    for pp = 1:length(params(mm).names)
        subplot(length(params(mm).names), 1, pp); hold on;
        for rr = 1:length(rois)
            pv = params(mm).values{rr, pp};
            p_avg = mean(pv); p_sem = std(pv)/sqrt(length(pv) - 1);
            p_min = p_avg - p_sem; p_max = p_avg + p_sem;
            p_xvec = [-.4 .4 .4 -.4] + rr;
            p_yvec = [p_min p_min p_max p_max];
            patch(p_xvec, p_yvec, rcols(rr, :), 'EdgeColor', ecols(rr, :));
            plot(p_xvec(1:2), [p_avg p_avg], 'Color', ecols(rr, :));
        end
        switch params(mm).names{pp}
            case 'tau1'
                ylims = [0 100];
            case 'tau2'
                ylims = [0 100];
            case 'epsilon'
                ylims = [0 .1];
            case 'sigma'
                ylims = [0 .1];
            case 'tau_s'
                ylims = [0 10];
            case 'tau_t'
                ylims = [0 10];
            case 'tau_d'
                ylims = [0 10];
        end
        axis tight; xlim([0 length(rois) + 1]); %ylim(ylims);
        title(params(mm).names{pp});
        set(gca, 'XTick', 1:length(rois), 'XTickLabel', roi_names, ...
            'TickDir', 'out', 'FontSize', 6);
    end
    fname = ['param_comparison_' models{mm} '_'];
    for ff = 1:size(fit_exps, 1)
        fname = [fname '_fit' [fit_exps{ff, :}]];
    end
    saveas(fig, fullfile(fig_dir, [fname '.fig']), 'fig');
end

end
