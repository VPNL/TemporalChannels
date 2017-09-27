function roi = pool_split_halves(rois, models, fit_exps, val_exps)

%% check inputs
if nargin < 1 || isempty(rois)
    rois = {'IOG-faces' 'pFus-faces' 'mFus-faces' ...
        'LOS-bodies' 'ITG-bodies' 'MTG-bodies' 'OTS-bodies' ...
        'IOS-characters' 'pOTS-characters' 'mOTS-characters'};
else
    rois = force_cell(strrep(rois, '_', '-'));
end

if nargin < 2 || isempty(models)
    models = {'glm' 'htd' 'balloon' ...
        'cts-pow' 'cts-div' 'dcts' ...
        '2ch' '2ch-rect' '2ch-sqrt-rect' ...
        '2ch-pow' '2ch-div' '2ch-dcts' '2ch-opt'};
else
    models = force_cell(models);
end

if nargin < 3 || isempty(fit_exps)
    fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
else
    fit_exps = force_cell(fit_exps);
end

if nargin < 4 || isempty(val_exps)
    val_exps = {'ExpAe' 'ExpBe' 'ExpCe'; 'ExpAo' 'ExpBo' 'ExpCo'};
else
    val_exps = force_cell(val_exps);
end

%% specify paths to files in results/figures directories
project_dir = fullfile(RAID, 'projects', 'CategoryChannels');
model_dir = fullfile(project_dir, 'model', 'TemporalChannels');
res_dir = fullfile(model_dir, 'results');
fig_dir = fullfile(model_dir, 'figures');
% generate stem of data filenames
fstems = cell(1, size(fit_exps, 1));
for ff = 1:size(fit_exps, 1)
    fstem = ['_fit' [fit_exps{ff, :}]];
    for ee = 1:size(val_exps, 2)
        fstem = [fstem '_val' val_exps{ff, ee}];
    end
    fstems{ff} = [fstem '.mat'];
end

%% load data for all ROIs, models, experiments, and data folds
for rr = 1:length(rois)
    for mm = 1:length(models)
        for ff = 1:size(fit_exps, 1)
            fname = [rois{rr} '_' models{mm} fstems{ff}];
            fold(ff) = load(fullfile(res_dir, fname));
        end
        roi = pool_across_folds(fold(1).roi, fold(2).roi);
        fname = [rois{rr} '_' models{mm} '_fitExpAExpBExpC_split_half.mat'];
        save(fullfile(res_dir, fname), 'roi', '-v7.3');
    end
end



end
