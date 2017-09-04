function [rois, models] = grid_search(roi_init, model_init, session_num, nseeds)
% Performs a grid search to find seeds for model optimization. 
% 
% INPUTS
%   1) roi_init: ROI object fitted with default model parameters
%   2) model_init: ModelTS object with default parameters
%   3) session_num: session number to perform grid search in
%   4) nseeds: number of model seed points to return
% 
% OUTPUTS
%   1) rois: updated ROI object fit with each of the top nseeds models
%   2) models: updated ModelTS object storing model for each seed point
% 
% AS 5/2017

session = roi_init.sessions{session_num};
param_names = fieldnames(model_init.params);
if length(param_names) > 2
    grid_size = 12;
else
    grid_size = 30;
end
params_lims = cell(length(param_names), 1);
params_grid = cell(length(param_names), 1);
for pp = 1:length(param_names)
    switch param_names{pp}
        case 'epsilon'
            params_lims{pp} = [.01 .99];
        case 'tau1'
            params_lims{pp} = [10 500];
        case 'tau2'
            params_lims{pp} = [10 500];
        case 'sigma'
            params_lims{pp} = [.01 .5];
    end
    params_grid{pp} = linspace(params_lims{pp}(1), params_lims{pp}(2), grid_size);
end

% generate all possible permutations of grid parameter values
nparams = length(param_names);
params_list = cell(nparams,1);
[params_list{:}] = ndgrid(params_grid{:});
params_list = reshape(cat(nparams+1,params_list{:}),[],nparams);
if sum(strcmp('tau1', param_names)) && sum(strcmp('tau2', param_names))
    plist = []; cnt = 0;
    tau1_idx = find(strcmp('tau1', param_names));
    tau2_idx = find(strcmp('tau2', param_names));
    for mm = 1:size(params_list, 1)
        if params_list(mm, tau1_idx) < params_list(mm, tau2_idx)
            cnt = cnt+1;
            plist(cnt, :) = params_list(mm, :);
        end
    end
    params_list = plist;
end
niters = size(params_list, 1);

% generate model grid for search only if necessary
fname = ['grid_search_models_' model_init.type '_fit' [roi_init.experiments{:}] '.mat'];
fpath = fullfile(session, 'Stimuli', fname);
ss_model = ModelTS(model_init.type, model_init.experiments, session);
ss_model = code_stim(ss_model);
if ~(exist(fpath, 'file') == 2)
    fprintf('Generating %s model grid for %s...\n', model_init.type, roi_init.session_ids{session_num});
    run_preds = {};
    for ii = 1:niters
        fprintf([num2str(ii) ' '])
        ii_model = ss_model;
        for pp = 1:length(param_names)
            ii_model.params.(param_names{pp}) = {params_list(ii, pp)};
            ii_model = update_param(ii_model, param_names{pp}, 0);
        end
        ii_model = pred_runs(ii_model);
        run_preds(:, ii) = ii_model.run_preds;
        if rem(ii, 20) == 0
            fprintf('\n')
        end
    end
    fprintf('\n')
    save(fpath, 'ii_model', 'run_preds', 'params_list', '-v7.3');
else
    load(fpath);
end

% fit models and find the best fitting parameter sets for each session
fprintf('Performing %s grid search for %s...\n', model_init.type, roi_init.session_ids{session_num});
ss_roi = ROI(roi_init.name, roi_init.experiments, session);
ss_roi = tc_runs(ss_roi);
ss_roi = tc_trials(ss_roi, ss_model);

% test
var_exp = {};
for ii = 1:niters
    fprintf([num2str(ii) ' '])
    ii_roi = ss_roi; ii_model = ss_model;
    for pp = 1:length(param_names)
        ii_model.params.(param_names{pp}) = {params_list(ii, pp)};
        ii_model = update_param(ii_model, param_names{pp}, 0);
    end
    ii_model.run_preds = run_preds(:, ii);
    ii_roi = tc_fit(ii_roi, ii_model);
    var_exp{ii} = [ii_roi.model.varexp{:}];
    if rem(ii, 20) == 0
        fprintf('\n')
    end
end
fprintf('\n')
[~, model_idxs] = sort(cell2mat(var_exp), 2, 'descend');


% search_rois = ROI(roi_init.name, roi_init.experiments, repmat({session}, 1, niters));
% search_rois = select_sessions(search_rois);
% search_rois.runs = repmat(ss_roi.runs, 1, niters);
% search_rois.baseline = repmat(ss_roi.baseline, 1, niters);
% search_rois = tc_fit(search_rois, search_models);
% [~, model_idxs] = sort([search_rois.model.varexp{:}], 2, 'descend');

% initialize output objects with correct number of seed points
for seed = 1:nseeds
    rois(seed) = ss_roi;
    models(seed) = ss_model;
    for pp = 1:length(param_names)
        % opt_param = search_models.params.(param_names{pp}){model_idxs(seed)};
        models(seed).params.(param_names{pp}){1} = params_list(seed, pp);
    end
end

% update IRFs and output results and parameters of best models
for seed = 1:nseeds
    for pp = 1:length(param_names)
        models(seed) = update_param(models(seed), param_names{pp}, 0);
    end
    models(seed) = pred_runs(models(seed));
    models(seed) = pred_trials(models(seed));
    rois(seed) = tc_fit(rois(seed), models(seed));
    rois(seed) = tc_pred(rois(seed), models(seed));
end

end
