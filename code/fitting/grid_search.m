function [roi, model] = grid_search(roi_init, model_init, nseeds)
% Performs a grid search to find starting seeds for model optimization. 
% 
% INPUTS
%   1) roi_init: ROI object fitted with default model parameters
%   2) model_init: ModelTS object with default parameters
%   3) nseeds: number of model seed points to return
% 
% OUTPUTS
%   1) roi: updated ROI object fit with each of the top nseeds models
%   2) model: updated ModelTS object storing model for each seed point
% 
% AS 5/2017

grid_size = 10;
sessions = roi_init.sessions; nsess = length(sessions);
param_names = fieldnames(model_init.params);
params_lims = cell(length(param_names), 1);
params_grid = cell(length(param_names), 1);
for pp = 1:length(param_names)
    switch param_names{pp}
        case 'e'
            params_lims{pp} = [.01 .99];
        case 'tau1'
            params_lims{pp} = [10 299];
        case 'tau2'
            params_lims{pp} = [11 300];
        case 'sigma'
            params_lims{pp} = [.01 .2];
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

% generate model grid for search only if necessary
sessions = roi_init.sessions; nsess = length(sessions);
fname = ['model_grid_' model_init.type '_' [roi_init.experiments{:}] '.mat'];
for ss = 1:length(sessions)
    [~, sess_id] = fileparts(sessions{ss});
    fpath = fullfile(roi_init.project_dir, 'data', sess_id, 'Stimuli', fname);
    if ~(exist(fpath, 'file') == 2)
        fprintf('Generating model grid for session %d of %d: ', ss, nsess);
        ss_model = ModelTS(model_init.type, model_init.experiments, sessions{ss});
        ss_model = code_stim(ss_model);
        for mm = 1:size(params_list, 1)
            fprintf('%d ', mm);
            models(mm) = ss_model;
            for pp = 1:length(param_names)
                models(mm).params.(param_names{pp}) = repmat({params_list(mm, pp)}, 1, 1);
                models(mm) = update_param(models(mm), param_names{pp}, 0);
            end
            models(mm) = run_preds(models(mm));
            models(mm) = trial_preds(models(mm));
        end
        fprintf('\n');
        save(fpath, 'models', '-v7.3');
        clear models;
    end
end

% initialize output objects with correct number of seed points
for seed = 1:nseeds
    roi(seed) = roi_init;
    model(seed) = model_init;
end

% fit models and find the best fitting parameter sets for each session
for ss = 1:length(sessions)
    fprintf('Performing grid search for session %d of %d... \n', ss, nsess);
    var_exp_grid = zeros(size(params_list, 1), 1);
    load(fullfile(sessions{ss}, 'Stimuli', fname));
    [~, session_id] = fileparts(sessions{ss});
    ss_roi = ROI(roi_init.name, roi_init.experiments, session_id);
    ss_roi = tc_runs(ss_roi);
    ss_roi = tc_trials(ss_roi, models(1));
    for mm = 1:size(params_list, 1)
        rois(mm) = ss_roi;
        rois(mm) = tc_fit(rois(mm), models(mm));
        rois(mm) = tc_pred(rois(mm), models(mm));
        var_exp_grid(mm) = [rois(mm).model.varexp{:}];
    end
    [~, model_idxs] = sort(var_exp_grid, 1, 'descend');
    for seed = 1:nseeds
        for pp = 1:length(param_names)
            max_param = models(model_idxs(seed)).params.(param_names{pp}){1};
            model(seed).params.(param_names{pp}){ss} = max_param;
        end
    end
    clear rois; clear models;
end

% update IRFs and output results and parameters of best nseeds models
for seed = 1:nseeds
    for pp = 1:length(param_names)
        model(seed) = update_param(model(seed), param_names{pp}, 0);
    end
    model(seed) = run_preds(model(seed));
    model(seed) = trial_preds(model(seed));
    roi(seed) = tc_fit(roi(seed), model(seed));
    roi(seed) = tc_pred(roi(seed), model(seed));
end


end
