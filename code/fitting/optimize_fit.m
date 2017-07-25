function [roi, model] = optimize_fit(roi_init, model_init)
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

nseeds = length(roi_init);
niterations = 20; step_decay = 0.9;
params_names = fieldnames(model_init(1).params);
sessions = model_init(1).sessions; nsess = length(sessions);
params_init = cell(length(params_names), nsess);
step_sizes_init = zeros(length(params_names), 1);
for pp = 1:length(params_names)
    params_init(pp, :) = model_init(1).params.(params_names{pp});
    switch params_names{pp}
        case 'e'
            step_sizes_init(pp) = .01;
        case 'tau1'
            step_sizes_init(pp) = 10;
        case 'tau2'
            step_sizes_init(pp) = 10;
        case 'sigma'
            step_sizes_init(pp) = .01;
    end
end

fprintf('Performing iterative gradient descent...\n');
for mm = 1:nseeds
    fprintf(['  Seed ' num2str(mm) ' of ' num2str(nseeds) ': ']);
    step_sizes = step_sizes_init;
    for ii = 1:niterations
        fprintf([num2str(ii) ' ']);
        % loop through all parameters in each iteration
        for pp = 1:length(params_names)
            % store variance explained for current position
            var_exp_init = [roi_init(mm).model.varexp{:}];
            param_init = model_init(mm).params.(params_names{pp});
            % take a step in positive direction and fit model
            model_pos = update_param(model_init(mm), params_names{pp}, step_sizes(pp));
            param_pos = model_pos.params.(params_names{pp});
            model_pos = run_preds(model_pos);
            model_pos = trial_preds(model_pos);
            roi_pos = tc_fit(roi_init(mm), model_pos);
            roi_pos = tc_pred(roi_pos, model_pos);
            var_exp_pos = [roi_pos.model.varexp{:}];
            % take a step in negative direction and fit model
            model_neg = update_param(model_init(mm), params_names{pp}, -step_sizes(pp));
            param_neg = model_neg.params.(params_names{pp});
            model_neg = run_preds(model_neg);
            model_neg = trial_preds(model_neg);
            roi_neg = tc_fit(roi_init(mm), model_neg);
            roi_neg = tc_pred(roi_neg, model_neg);
            var_exp_neg = [roi_neg.model.varexp{:}];
            % find the best performing version of model for each session
            var_exp_ii = [var_exp_init; var_exp_pos; var_exp_neg];
            param_ii = [param_init; param_pos; param_neg];
            [~, opt_idxs] = max(var_exp_ii);
            for ss = 1:nsess
                param_new{ss} = param_ii{opt_idxs(ss), ss};
            end
            % update model_init to have best fitting parameters
            model_init(mm).params.(params_names{pp}) = param_new;
            model_init(mm) = update_param(model_init(mm), params_names{pp}, 0);
            model_init(mm) = run_preds(model_init(mm));
            model_init(mm) = trial_preds(model_init(mm));
            roi_init(mm) = tc_fit(roi_init(mm), model_init(mm));
            roi_init(mm) = tc_pred(roi_init(mm), model_init(mm));
        end
        step_sizes = step_sizes * step_decay;
    end
    fprintf('\n');
end

% store best fitting model
varexp = zeros(1, nseeds);
for mm = 1:nseeds
    varexp(mm) = mean([roi_init(mm).model.varexp{:}]);
end
[~, model_idx] = max(varexp);
roi = roi_init(mm); model = model_init(mm);

end
