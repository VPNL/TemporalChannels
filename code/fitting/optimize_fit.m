function [rois, models] = optimize_fit(roi_init, model_init)
% Performs an iterative gradient descent procedure to optimize model fit.
% 
% INPUTS
%   1) roi_init: ROI object fitted with seed points from grid seach
%   2) model_init: ModelTS object with models used to fit roi_init
% 
% OUTPUTS
%   1) roi: updated ROI object fit with optimized model parameters
%   2) model: updated ModelTS object used to fit roi
% 
% AS 5/2017

roi_iter = roi_init; model_iter = model_init;
nseeds = length(roi_iter);
niter = 20; step_decay = 0.9;
param_names = fieldnames(model_iter(1).params);
sessions = model_iter(1).sessions; nsess = length(sessions);
params_init = cell(length(param_names), nsess);
step_sizes_init = zeros(length(param_names), 1);
for pp = 1:length(param_names)
    params_init(pp, :) = model_iter(1).params.(param_names{pp});
    switch param_names{pp}
        case 'epsilon'
            step_sizes_init(pp) = .01;
        case 'tau1'
            step_sizes_init(pp) = 10;
        case 'tau2'
            step_sizes_init(pp) = 10;
        case 'sigma'
            step_sizes_init(pp) = .01;
        case 'tau_s'
            step_sizes_init(pp) = .5;
        case 'tau_t'
            step_sizes_init(pp) = .5;
        case 'tau_d'
            step_sizes_init(pp) = .5;
    end
end

fprintf('Performing iterative parameter updates for %s...\n', roi_init(1).session_ids{1});
for mm = 1:nseeds
    fprintf('  Seed %d of %d: ', mm, nseeds);
    step_sizes = step_sizes_init;
    for ii = 1:niter
        fprintf([num2str(ii) ' ']);
        % loop through all parameters in each iteration
        for pp = 1:length(param_names)
            % store variance explained for current position
            var_exp_init = [roi_iter(mm).model.varexp{:}];
            param_init = model_iter(mm).params.(param_names{pp});
            % take a step in positive direction and fit model
            model_pos = update_param(model_iter(mm), param_names{pp}, step_sizes(pp));
            param_pos = model_pos.params.(param_names{pp});
            model_pos = pred_runs(model_pos);
            model_pos = pred_trials(model_pos);
            roi_pos = tc_fit(roi_iter(mm), model_pos);
            var_exp_pos = [roi_pos.model.varexp{:}];
            % take a step in negative direction and fit model
            model_neg = update_param(model_iter(mm), param_names{pp}, -step_sizes(pp));
            param_neg = model_neg.params.(param_names{pp});
            model_neg = pred_runs(model_neg);
            model_neg = pred_trials(model_neg);
            roi_neg = tc_fit(roi_iter(mm), model_neg);
            var_exp_neg = [roi_neg.model.varexp{:}];
            % find the best performing version of model for each session
            var_exp_ii = [var_exp_init var_exp_pos var_exp_neg];
            param_ii = [param_init param_pos param_neg];
            [~, opt_idxs] = max(var_exp_ii);
            param_new{1} = param_ii{opt_idxs(1)};
            % update model_init to have best fitting parameters
            model_iter(mm).params.(param_names{pp}) = param_new;
            model_iter(mm) = update_param(model_iter(mm), param_names{pp}, 0);
            model_iter(mm) = pred_runs(model_iter(mm));
            roi_iter(mm) = tc_fit(roi_iter(mm), model_iter(mm));
        end
        step_sizes = step_sizes * step_decay;
    end
    fprintf('\n');
end

% store best fitting models
varexp = zeros(1, nseeds);
for mm = 1:nseeds
    varexp(mm) = roi_iter(mm).model.varexp{:};
end
[~, model_idx] = max(varexp);
rois = roi_iter(model_idx);
models = model_iter(model_idx);
models = pred_trials(models);


end
