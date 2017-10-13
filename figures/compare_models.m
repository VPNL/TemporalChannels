function compare_models(rois, models)
% Generates model comparison figure for all ROIs by plotting R^2 in each
% experiment in val_exps for a model fit across all fit_exps data.
% When multiple rows of fit_exps and val_exps are passed, performance is
% averaged across both splits of the data.
%
% Inputs
%   1) rois: cell array of ROI names
%   2) models: cell array of model names
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
    models = {'1ch-lin' '2ch-lin-htd' '1ch-balloon' ...
        '1ch-pow' '1ch-div' '1ch-dcts' '1ch-exp' '1ch-cexp' ...
        '2ch-lin-quad' '2ch-lin-rect' ...
        '2ch-pow-quad' '2ch-pow-rect' ...
        '2ch-exp-quad' '2ch-exp-rect' ...
        '2ch-cexp-quad' '2ch-cexp-rect' ...
        '3ch-lin-quad-exp' '3ch-lin-rect-exp' ...
        '3ch-exp-quad-exp' '3ch-exp-rect-exp' ...
        '3ch-cexp-quad-exp' '3ch-cexp-rect-exp'};
    mnames = {'GLM' 'HTD' 'Balloon' ...
        'CTS-p' 'CTS-n' 'dCTS' 'Adapt' 'cAdapt' ...
        '2ch-lq' '2ch-lr' ...
        '2ch-pq' '2ch-pr' ...
        '2ch-eq' '2ch-er' ...
        '2ch-ceq' '2ch-cer' ...
        '3ch-lqe' '3ch-lre' ...
        '3ch-eqe' '3ch-ere' ...
        '3ch-ceqe' '3ch-cere'};
    mcols = [repmat([1 1 1], 3, 1); ...
        repmat([.8 .8 .8], 5, 1); ...
        repmat([.6 .6 .6], 8, 1); ...
        repmat([.4 .4 .4], 6, 1)];
else
    [models, mnames] = deal(force_cell(models));
    mcols = repmat([0 0 0], length(models), 1);
end

%% specify paths to files in results/figures directories
project_dir = fullfile(RAID, 'projects', 'CategoryChannels');
model_dir = fullfile(project_dir, 'model', 'TemporalChannels');
res_dir = fullfile(model_dir, 'results');
fig_dir = fullfile(model_dir, 'figures', 'model_comparison');

%% load data for all ROIs, models, experiments, and data folds
var_exp = {}; noise_ceils = {};
for rr = 1:length(rois)
    for mm = 1:length(models)
        fname = [rois{rr} '_' models{mm} '_split-half.mat'];
        load(fullfile(res_dir, fname));
        for ee = 1:length(roi) - 1
            var_exp{rr, mm, ee} = [roi(ee + 1).model.varexp{:}];
            if mm == 1
                noise_ceils{rr, ee} = [roi(ee + 1).noise_ceils{:}];
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
fig = tch_fig('Model comparison', [.05 .05 .8 .8]);
for ee = 1:length(val_names)
    subplot(length(val_names), 1, ee); hold on; tch_set_axes(gca);
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
        ve = cell2mat(var_exp(rr, :, ee)')';
        tch_plot_bar(ve, mcols, xlims(rr, 1));
        text(mean(xlims(rr, :)), ylims(2), rois{rr}, 'FontSize', 6, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
        if rr == 1 && ee == 1
            for mm = 1:length(models)
                text(mm + xlims(rr, 1), .35, mnames{mm}, 'FontSize', 6, 'Rotation', 90, ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
            end
        end
    end
end
fname = ['model_comparison_fit' [fit_exps{1, :}] '_val' [val_exps{1, :}]];
saveas(fig, fullfile(fig_dir, [fname '.fig']), 'fig');
print(fullfile(fig_dir, [fname '.eps']), '-depsc2', '-painters');
close all;

end
