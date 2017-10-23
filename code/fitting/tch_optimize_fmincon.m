function [roi, model] = tch_optimize_fmincon(roi, model, fit_exps)
% Nonlinear optimization prodcedure using built-in fmincon algorithm.

param_names = fieldnames(model.params); sessions = roi.sessions;
for ss = 1:length(sessions)
    fname_opt = ['optimization_results_' model.type '_fit' [fit_exps{:}] '.mat'];
    fpath_opt = fullfile(sessions{ss}, 'ROIs', roi.name, fname_opt);
    % load optimization results if saved, otherwise compute
    if exist(fpath_opt, 'file') == 2
        load(fpath_opt); fprintf('Loading gradient descent results. \n');
        omodel = code_stim(tchModel(model.type, roi.experiments, sessions{ss}));
        omodel.normT = model.normT; omodel.normD = model.normD;
        for pp = 1:length(param_names)
            pn = param_names{pp}; omodel.params.(pn){1} = params.(pn){1};
        end
        omodel = update_param(omodel, pn, 0);
    else
        sroi = tch_runs(tchROI(roi.name, roi.experiments, sessions{ss}));
        omodel = code_stim(tchModel(model.type, roi.experiments, sessions{ss}));
        omodel = pred_runs(omodel); omodel = pred_trials(omodel);
        sroi = tch_trials(sroi, omodel); sroi = tch_fit(sroi, omodel);
        npreds = size(sroi.model.betas{1}, 2);
        fmin_options = optimoptions('fmincon', 'Display', 'iter', ...
            'StepTolerance', 1e-2, 'UseParallel', true);
        switch model.type
            case '1ch-pow'
                obj_fun = tch_obj_fun_1ch_pow(sroi, omodel);
                x_init = [50 .1 sroi.model.betas{1}];
                lb = [10 .01 -Inf(1, npreds)];
                ub = [1000 1 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau1{1} = x_opt(1);
                params.epsilon{1} = x_opt(2);
            case '1ch-exp'
                obj_fun = tch_obj_fun_1ch_exp(sroi, omodel);
                x_init = [10000 sroi.model.betas{1}];
                lb = [100 -Inf(1, npreds)];
                ub = [60000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_ae{1} = x_opt(1);
            case '2ch-pow-quad'
                obj_fun = tch_obj_fun_2ch_pow_quad(sroi, omodel);
                x_init = [.1 sroi.model.betas{1}];
                lb = [.01 -Inf(1, npreds)];
                ub = [1 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.epsilon{1} = x_opt(1);
            case '2ch-pow-rect'
                obj_fun = tch_obj_fun_2ch_pow_rect(sroi, omodel);
                x_init = [.1 sroi.model.betas{1}];
                lb = [.01 -Inf(1, npreds)];
                ub = [1 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.epsilon{1} = x_opt(1);
            case '2ch-exp-quad'
                obj_fun = tch_obj_fun_2ch_exp_quad(sroi, omodel);
                x_init = [10000 sroi.model.betas{1}];
                lb = [100 -Inf(1, npreds)];
                ub = [60000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_ae{1} = x_opt(1);
            case '2ch-exp-rect'
                obj_fun = tch_obj_fun_2ch_exp_rect(sroi, omodel);
                x_init = [10000 sroi.model.betas{1}];
                lb = [100 -Inf(1, npreds)];
                ub = [60000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_ae{1} = x_opt(1);
            case '3ch-lin-quad-exp'
                obj_fun = tch_obj_fun_3ch_lin_quad_exp(sroi, omodel);
                x_init = [1000 sroi.model.betas{1}];
                lb = [10 -Inf(1, npreds)];
                ub = [12000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_de{1} = x_opt(1);
            case '3ch-lin-rect-exp'
                obj_fun = tch_obj_fun_3ch_lin_rect_exp(sroi, omodel);
                x_init = [1000 sroi.model.betas{1}];
                lb = [10 -Inf(1, npreds)];
                ub = [12000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_de{1} = x_opt(1);
            case '3ch-pow-quad-exp'
                obj_fun = tch_obj_fun_3ch_pow_quad_exp(sroi, omodel);
                x_init = [.1 1000 sroi.model.betas{1}];
                lb = [.01 10 -Inf(1, npreds)];
                ub = [1 12000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.epsilon{1} = x_opt(1);
                params.tau_de{1} = x_opt(2);
            case '3ch-pow-rect-exp'
                obj_fun = tch_obj_fun_3ch_pow_rect_exp(sroi, omodel);
                x_init = [.1 1000 sroi.model.betas{1}];
                lb = [.01 10 -Inf(1, npreds)];
                ub = [1 12000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.epsilon{1} = x_opt(1);
                params.tau_de{1} = x_opt(2);
            case '3ch-exp-quad-exp'
                obj_fun = tch_obj_fun_3ch_exp_quad_exp(sroi, omodel);
                x_init = [10000 10000 sroi.model.betas{1}];
                lb = [100 10 -Inf(1, npreds)];
                ub = [60000 12000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_ae{1} = x_opt(1);
                params.tau_de{1} = x_opt(2);
            case '3ch-exp-rect-exp'
                obj_fun = tch_obj_fun_3ch_exp_rect_exp(sroi, omodel);
                x_init = [10000 10000 sroi.model.betas{1}];
                lb = [100 10 -Inf(1, npreds)];
                ub = [60000 12000 Inf(1, npreds)];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], lb, ub, [], fmin_options);
                params.tau_ae{1} = x_opt(1);
                params.tau_de{1} = x_opt(2);
        end
        for pp = 1:length(param_names)
            pn = param_names{pp};
            omodel.params.(pn){1} = params.(pn){1};
        end
        omodel = update_param(omodel, param_names{pp}, 0);
        save(fpath_opt, 'params', '-v7.3');
    end
    % copy optimized parameters from session to group model objects
    for pp = 1:length(param_names)
        pn = param_names{pp}; model.params.(pn){ss} = omodel.params.(pn){1};
    end
end
model = update_param(model, param_names{pp}, 0); model = norm_model(model, 1);
model = pred_runs(model); model = pred_trials(model);
[roi, model] = tch_fit(roi, model, 0);

end
