function [oroi, omodel] = tch_optimize_fit(iroi, imodel)
% Performs a gradient descent procedure with iterative coordinate updates
% to optimize model fit for a given ROI.
% 
% INPUTS
%   1) iroi: initial ROI object fitted with seed points from grid seach
%   2) imodel: inital tchModel object with model used to fit iroi
% 
% OUTPUTS
%   1) oroi: updated tchROI object fit with optimized model parameters
%   2) omodel: updated tchModel object with parameters used to fit oroi
% 
% AS 5/2017

% setup variables to store fit iterations and track steps
roi_iter = iroi; model_iter = imodel; 
nseeds = length(roi_iter); niter = 20; step_decay = 0.9;
param_names = fieldnames(model_iter(1).params);
sessions = model_iter(1).sessions; nsess = length(sessions);
params_init = cell(length(param_names), nsess);
step_sizes_init = zeros(length(param_names), 1);

% get inital step size for each parameter
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
        case 'kappa'
            step_sizes_init(pp) = .1;
        case 'tau_ae'
            step_sizes_init(pp) = 500;
    end
end

fprintf('Performing coordinate updates for %s...\n', iroi(1).session_ids{1});
for mm = 1:nseeds
    step_sizes = step_sizes_init; fprintf('  Seed %d of %d: ', mm, nseeds);
    for ii = 1:niter
        fprintf([num2str(ii) ' ']);
        % loop through all parameters in each iteration
        for pp = 1:length(param_names)
            % store variance explained for current parameter coordinates
            var_exp_init = [roi_iter(mm).model.varexp{:}];
            param_init = model_iter(mm).params.(param_names{pp});
            pn = param_names{pp};
            % take a step in positive direction and fit model
            model_pos = tch_update_param(model_iter(mm), pn, step_sizes(pp));
            param_pos = model_pos.params.(param_names{pp});
            model_pos = pred_runs(model_pos);
            model_pos = pred_trials(model_pos);
            roi_pos = tch_fit(roi_iter(mm), model_pos);
            var_exp_pos = [roi_pos.model.varexp{:}];
            % take a step in negative direction and fit model
            model_neg = tch_update_param(model_iter(mm), pn, -step_sizes(pp));
            param_neg = model_neg.params.(param_names{pp});
            model_neg = pred_runs(model_neg);
            model_neg = pred_trials(model_neg);
            roi_neg = tch_fit(roi_iter(mm), model_neg);
            var_exp_neg = [roi_neg.model.varexp{:}];
            % find the best performing parameter coordinates
            var_exp_ii = [var_exp_init var_exp_pos var_exp_neg];
            param_ii = [param_init param_pos param_neg];
            [~, opt_idxs] = max(var_exp_ii);
            param_new{1} = param_ii{opt_idxs(1)};
            % update model_init to have best fitting parameters
            model_iter(mm).params.(param_names{pp}) = param_new;
            model_iter(mm) = tch_update_param(model_iter(mm), param_names{pp}, 0);
            model_iter(mm) = pred_runs(model_iter(mm));
            roi_iter(mm) = tch_fit(roi_iter(mm), model_iter(mm));
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
[~, opt_idx] = max(varexp); oroi = roi_iter(opt_idx);
omodel = model_iter(opt_idx); omodel = pred_trials(omodel);

end
