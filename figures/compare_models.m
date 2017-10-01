function compare_models(rois, models, fit_exps, val_exps)
% Generates model comparison figure for all ROIs by plotting R^2 in each
% experiment in val_exps for a model fit across all fit_exps data.
% When multiple rows of fit_exps and val_exps are passed, performance is
% averaged across both splits of the data.
%
% Inputs
%   1) rois: cell array of ROI names (e.g., {'IOG-faces' 'pFus-faces'})
%   2) models: cell array of model names (e.g., {'standard' '2ch-dcts'})
%   3) fit_exps: cell array of experiments used to fit the model
%   4) val_exps: cell array of experiments to validate the model
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

if nargin < 2 || isempty(models)
    models = {'glm' 'htd' 'balloon' ...
        'cts-pow' 'cts-div' 'dcts' ...
        '2ch' '2ch-rect' '2ch-pow' ...
        '2ch-opt' '3ch' '3ch-opt'};
    mnames = {'GLM' 'HTD' 'Balloon' ...
        'CTS-p' 'CTS-n' 'dCTS' ...
        '2ch' '2ch-r' '2ch-p' ...
        '2ch-o' '3ch' '3ch-o'};
    mcols = [1 1 1; 1 1 1; 1 1 1; ...
        .75 .75 .75; .75 .75 .75; .75 .75 .75; ...
        .5 .5 .5; .5 .5 .5; .5 .5 .5; ...
        .25 .25 .25; .25 .25 .25; .25 .25 .25];
else
    models = force_cell(models);
    mnames = models;
    mcols = repmat([0 0 0], length(models), 1);
end

if nargin < 3 || isempty(fit_exps)
    fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
else
    fit_exps = force_cell(fit_exps);
end

if nargin < 4 || isempty(val_exps)
    val_exps = {'ExpAe' 'ExpBe' 'ExpCe'; 'ExpAo' 'ExpBo' 'ExpCo'};
    val_names = {'Exp1' 'Exp2' 'Exp3'};
else
    val_exps = force_cell(val_exps);
    val_names = val_exps;
end

%% specify paths to files in results/figures directories
project_dir = fullfile(RAID, 'projects', 'CategoryChannels');
model_dir = fullfile(project_dir, 'model', 'TemporalChannels');
res_dir = fullfile(model_dir, 'results', 'pooled');
fig_dir = fullfile(model_dir, 'figures', 'model_comparison');

%% load data for all ROIs, models, experiments, and data folds
var_exp = {}; noise_ceils = {};
for rr = 1:length(rois)
    for mm = 1:length(models)
        fname = [rois{rr} '_' models{mm} '_split_half.mat'];
        if exist(fullfile(res_dir, fname), 'file') == 2
            load(fullfile(res_dir, fname));
            for ee = 1:size(val_exps, 2)
                var_exp{rr, mm, ee} = [roi(ee + 1).model.varexp{:}];
                if mm == 1
                    noise_ceils{rr, ee} = [roi(ee + 1).noise_ceils{:}];
                end
            end
        else
            for ee = 1:size(val_exps, 2)
                var_exp{rr, mm, ee} = zeros(size(var_exp{rr, 1, ee}));
            end
        end
    end
end

%% plot model comparison for each experiment and ROI averaged across folds
gap = 1; xcnt = 1; xlims = zeros(length(val_names), 2); ylims = [-.4 1.001];
for rr = 1:length(rois)
    xlims(rr, 1) = xcnt; xlims(rr, 2) = xcnt + length(models) + 1;
    xcnt = xcnt + length(models) + 1 + gap;
end
fig = figTS('Model comparison', [.05 .05 .9 .9]);
for ee = 1:length(val_names)
    subplot(length(val_names), 1, ee); hold on;
    tc_axes(val_names{ee}, '', 'Cross-validated R^2', xlims, ylims);
    for rr = 1:length(rois)
        nc = noise_ceils{rr, ee};
        nc_avg = mean(nc); nc_sem = std(nc) / sqrt(length(nc) - 1);
        nc_min = nc_avg - nc_sem; nc_max = nc_avg + nc_sem;
        nc_xvec = [.5 length(models) + .5] + xlims(rr, 1);
        nc_xvec = [nc_xvec fliplr(nc_xvec)];
        nc_yvec = [nc_min nc_min nc_max nc_max];
        patch(nc_xvec, nc_yvec, [.9 .9 .9], 'EdgeColor', [.8 .8 .8]);
        plot(nc_xvec(1:2), [nc_avg nc_avg], 'k--');
        ve = cell2mat(var_exp(rr, :, ee)')'; barTS(ve, mcols, xlims(rr, 1));
        text(mean(xlims(rr, :)), ylims(2), rois{rr}, 'FontSize', 7, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        if rr == 1 && ee == 1
            for mm = 1:length(models)
                text(mm + xlims(rr, 1), .35, mnames{mm}, 'FontSize', 6, ...
                    'HorizontalAlignment', 'left', ...
                    'VerticalAlignment', 'middle', ...
                    'Rotation', 90);
            end
        end
    end
end
fname = ['model_comparison_fit' [fit_exps{1, :}] '_val' [val_exps{1, :}]];
saveas(fig, fullfile(fig_dir, [fname '.fig']), 'fig');
print(fullfile(fig_dir, [fname '_early.eps']), '-depsc2', '-painters');
close all;

end
