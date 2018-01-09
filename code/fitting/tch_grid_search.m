function [rois, models] = tch_grid_search(iroi, imodel, session_num, nseeds)
% Performs a grid search to find seeds for further model optimization. 
% 
% INPUTS
%   1) iroi: tchROI object fitted with default model parameters
%   2) imodel: tchModel object with default parameters
%   3) session_num: session number to perform grid search in
%   4) nseeds: number of model seed points to return
% 
% OUTPUTS
%   1) rois: updated tchROI object fit with each of the top nseeds models
%   2) models: updated tchModel object storing model for each seed point
% 
% AS 5/2017

% set grid size for search
session = iroi.sessions{session_num};
param_names = fieldnames(imodel.params);
if length(param_names) > 2
    grid_size = 12;
else
    grid_size = 30;
end
[plims, params_grid] = deal(cell(length(param_names), 1));

% specify min and max values to search for each parameter
for pp = 1:length(param_names)
    switch param_names{pp}
        case 'epsilon'
            plims{pp} = [.01 .99];
        case 'tau1'
            plims{pp} = [10 500];
        case 'tau2'
            plims{pp} = [10 500];
        case 'sigma'
            plims{pp} = [.01 .5];
        case 'tau_s'
            plims{pp} = [1 50];
        case 'tau_t'
            plims{pp} = [1 50];
        case 'tau_d'
            plims{pp} = [1 50];
        case 'kappa'
            plims{pp} = [1 10];
        case 'tau_ae'
            plims{pp} = [100 60000];
        case 'tau_pe'
            plims{pp} = [10 12000];
    end
    params_grid{pp} = linspace(plims{pp}(1), plims{pp}(2), grid_size);
end

% generate all possible permutations of grid parameter values
nparams = length(param_names); params_list = cell(nparams, 1);
[params_list{:}] = ndgrid(params_grid{:});
params_list = reshape(cat(nparams + 1, params_list{:}), [], nparams);

% remove invalid models (when tau2 < tau1)
if sum(strcmp('tau1', param_names)) && sum(strcmp('tau2', param_names))
    plist = []; cnt = 0;
    tau1_idx = find(strcmp('tau1', param_names));
    tau2_idx = find(strcmp('tau2', param_names));
    for mm = 1:size(params_list, 1)
        if params_list(mm, tau1_idx) < params_list(mm, tau2_idx)
            cnt = cnt+1; plist(cnt, :) = params_list(mm, :);
        end
    end
    params_list = plist;
end
niters = size(params_list, 1);

% generate model grid only if necessary
fname = ['grid_search_models_' imodel.type '_fit' [iroi.experiments{:}] '.mat'];
fpath = fullfile(session, 'Stimuli', fname);
ss_model = tchModel(imodel.type, imodel.experiments, session);
ss_model.tr = imodel.tr; ss_model = code_stim(ss_model);
if ~(exist(fpath, 'file') == 2)
    fprintf('Generating %s model grid for %s...\n', ...
        imodel.type, iroi.session_ids{session_num});
    run_preds = {};
    for ii = 1:niters
        ii_model = ss_model; fprintf([num2str(ii) ' ']);
        for pp = 1:length(param_names)
            ii_model.params.(param_names{pp}) = {params_list(ii, pp)};
            ii_model = tch_update_param(ii_model, param_names{pp}, 0);
        end
        ii_model = pred_runs(ii_model);
        run_preds(:, ii) = ii_model.run_preds;
        if rem(ii, 20) == 0; fprintf('\n'); end
    end
    fprintf('\n')
    save(fpath, 'run_preds', 'params_list', '-v7.3');
else
    load(fpath);
end

% fit each model in grid and find the best fitting parameter sets
fprintf('Performing %s grid search for %s...\n', ...
    imodel.type, iroi.session_ids{session_num}); var_exp = {};
ss_roi = tchROI(iroi.name, iroi.experiments, session); ss_roi.tr = iroi.tr;
ss_roi = tch_runs(ss_roi); ss_roi = tch_trials(ss_roi, ss_model);
for ii = 1:niters
    ii_roi = ss_roi; ii_model = ss_model; fprintf([num2str(ii) ' ']);
    for pp = 1:length(param_names)
        ii_model.params.(param_names{pp}) = {params_list(ii, pp)};
        ii_model = tch_update_param(ii_model, param_names{pp}, 0);
    end
    ii_model.run_preds = run_preds(:, ii);
    ii_roi = tch_fit(ii_roi, ii_model);
    var_exp{ii} = [ii_roi.model.varexp{:}];
    if rem(ii, 20) == 0; fprintf('\n'); end
end
[~, model_idxs] = sort(cell2mat(var_exp), 2, 'descend'); fprintf('\n');

% initialize output objects with correct number of seed points
for seed = 1:nseeds
    rois(seed) = ss_roi; models(seed) = ss_model;
    for pp = 1:length(param_names)
        pidx = model_idxs(seed);
        models(seed).params.(param_names{pp}){1} = params_list(pidx, pp);
    end
end

% update IRFs and output results of nseeds best models
for seed = 1:nseeds
    for pp = 1:length(param_names)
        models(seed) = tch_update_param(models(seed), param_names{pp}, 0);
    end
    models(seed) = pred_runs(models(seed));
    models(seed) = pred_trials(models(seed));
    rois(seed) = tch_fit(rois(seed), models(seed));
    rois(seed) = tch_pred(rois(seed), models(seed));
end

end
