function [roi, model] = tch_optimize_fmincon(roi, model, fit_exps)
% Nonlinear optimization prodcedure using built-in fmincon algorithm.

param_names = fieldnames(model.params); sessions = roi.sessions;
if length(param_names) <= 3; fprec = 1e-6; end
if length(param_names) > 3; fprec = 1e-4; end
for ss = 1:length(sessions)
    fname_opt = ['optimization_results_' model.type '_fit' [fit_exps{:}] '.mat'];
    fpath_opt = fullfile(sessions{ss}, 'ROIs', roi.name, fname_opt);
    % load optimization results if saved, otherwise compute
    if exist(fpath_opt, 'file') == 2
        load(fpath_opt); fprintf('Loading gradient descent results. \n');
        omodel = code_stim(tchModel(model.type, roi.experiments, sessions{ss}));
        for pp = 1:length(param_names)
            pn = param_names{pp}; omodel.params.(pn){1} = params.(pn){1};
        end
        omodel = update_param(omodel, pn, 0);
    else
        sroi = tch_runs(tchROI(roi.name, roi.experiments, sessions{ss}));
        omodel = code_stim(tchModel(model.type, roi.experiments, sessions{ss}));
        omodel.normT = model.normT; omodel.normD = model.normD;
        omodel = pred_runs(omodel); omodel = pred_trials(omodel);
        sroi = tch_trials(sroi, omodel); sroi = tch_fit(sroi, omodel);
        pout = cell(1, length(param_names));
        fmin_options = optimoptions('fmincon', 'Display', 'off', ...
            'StepTolerance', fprec, 'UseParallel', true);
        fprintf('Optimizing parameters for %s ...\n', roi.session_ids{ss});
        switch model.type
            case '1ch-pow'
                obj_fun = tch_obj_fun_1ch_pow(sroi, omodel);
                x_init = [omodel.params.tau1{1} / 1000 omodel.params.epsilon{1}];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [.01 .001], [.5 1], [], fmin_options);
                params.tau1{1} = x_opt(1) * 1000;
                params.epsilon{1} = x_opt(2);
            case '1ch-div'
                obj_fun = tch_obj_fun_1ch_div(sroi, omodel);
                x_init = [omodel.params.tau1{1} / 1000 omodel.params.sigma{1}];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [.01 .001], [.5 1], [], fmin_options);
                params.tau1{1} = x_opt(1) * 1000;
                params.sigma{1} = x_opt(2);
            case '1ch-dcts'
                obj_fun = tch_obj_fun_1ch_dcts(sroi, omodel);
                x_init = [omodel.params.tau1{1} / 1000 ...
                    omodel.params.sigma{1} omodel.params.tau2{1} / 1000];
                x_opt = fmincon(obj_fun, x_init, [1 0 -1], 0, ...
                    [], [], [.01 .001 .01], [.5 1 1], [], fmin_options);
                params.tau1{1} = x_opt(1) * 1000;
                params.sigma{1} = x_opt(2);
                params.tau2{1} = x_opt(3) * 1000;
            case '1ch-exp'
                obj_fun = tch_obj_fun_1ch_exp(sroi, omodel);
                x_init = omodel.params.tau_ae{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 1, 60, [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
            case '2ch-pow-quad'
                obj_fun = tch_obj_fun_2ch_pow_quad(sroi, omodel);
                x_init = omodel.params.epsilon{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], .001, 1, [], fmin_options);
                params.epsilon{1} = x_opt(1);
            case '2ch-pow-rect'
                obj_fun = tch_obj_fun_2ch_pow_rect(sroi, omodel);
                x_init = omodel.params.epsilon{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], .001, 1, [], fmin_options);
                params.epsilon{1} = x_opt(1);
            case '2ch-exp-quad'
                obj_fun = tch_obj_fun_2ch_exp_quad(sroi, omodel);
                x_init = omodel.params.tau_ae{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 1, 60, [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
            case '2ch-exp-rect'
                obj_fun = tch_obj_fun_2ch_exp_rect(sroi, omodel);
                x_init = omodel.params.tau_ae{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 1, 60, [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
            case '3ch-lin-quad-exp'
                obj_fun = tch_obj_fun_3ch_lin_quad_exp(sroi, omodel);
                x_init = omodel.params.tau_de{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], .1, 12, [], fmin_options);
                params.tau_de{1} = x_opt(1) * 1000;
            case '3ch-lin-rect-exp'
                obj_fun = tch_obj_fun_3ch_lin_rect_exp(sroi, omodel);
                x_init = omodel.params.tau_de{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], .1, 12, [], fmin_options);
                params.tau_de{1} = x_opt(1) * 1000;
            case '3ch-pow-quad-exp'
                obj_fun = tch_obj_fun_3ch_pow_quad_exp(sroi, omodel);
                x_init = [omodel.params.epsilon{1} omodel.params.tau_de{1} / 1000];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [.001 .1], [1 12], [], fmin_options);
                params.epsilon{1} = x_opt(1);
                params.tau_de{1} = x_opt(2) * 1000;
            case '3ch-pow-rect-exp'
                obj_fun = tch_obj_fun_3ch_pow_rect_exp(sroi, omodel);
                x_init = [omodel.params.epsilon{1} omodel.params.tau_de{1} / 1000];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [.001 .1], [1 12], [], fmin_options);
                params.epsilon{1} = x_opt(1);
                params.tau_de{1} = x_opt(2) * 1000;
             case '3ch-exp-quad-exp'
                obj_fun = tch_obj_fun_3ch_exp_quad_exp(sroi, omodel);
                x_init = [omodel.params.tau_ae{1} omodel.params.tau_de{1}] / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [1 .1], [60 12], [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
                params.tau_de{1} = x_opt(2) * 1000;
            case '3ch-exp-rect-exp'
                obj_fun = tch_obj_fun_3ch_exp_rect_exp(sroi, omodel);
                x_init = [omodel.params.tau_ae{1} omodel.params.tau_de{1}] / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [1 .1], [60 12], [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
                params.tau_de{1} = x_opt(2) * 1000;
            case '2ch-lin-quad-opt'
                obj_fun = tch_obj_fun_2ch_lin_quad_opt(sroi, omodel);
                x_init = [omodel.params.tau_s{1} omodel.params.tau_t{1} omodel.params.kappa{1}]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 4 1], [50 50 10], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_t{1} = x_opt(2);
                params.kappa{1} = x_opt(3);
            case '2ch-lin-rect-opt'
                obj_fun = tch_obj_fun_2ch_lin_rect_opt(sroi, omodel);
                x_init = [omodel.params.tau_s{1} omodel.params.tau_t{1} omodel.params.kappa{1}]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 4 1], [50 50 10], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_t{1} = x_opt(2);
                params.kappa{1} = x_opt(3);
            case '3ch-lin-quad-exp-opt'
                obj_fun = tch_obj_fun_3ch_lin_quad_exp_opt(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_t = omodel.params.tau_t{1};
                tau_de = omodel.params.tau_de{1} / 1000;
                x_init = [tau_s tau_t omodel.params.kappa{1} tau_de]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 4 1 .1], [50 50 10 12], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_t{1} = x_opt(2);
                params.kappa{1} = x_opt(3);
                params.tau_de{1} = x_opt(4) * 1000;
            case '3ch-lin-rect-exp-opt'
                obj_fun = tch_obj_fun_3ch_lin_rect_exp_opt(sroi, omodel);
                tau_de = omodel.params.tau_de{1} / 1000;
                x_init = [omodel.params.tau_s{1} omodel.params.kappa{1} tau_de]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1 .1], [50 10 12], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.kappa{1} = x_opt(2);
                params.tau_de{1} = x_opt(3) * 1000;
        end
        for pp = 1:length(param_names)
            pn = param_names{pp}; pv = params.(pn){1}; pstr = num2str(pv, 3);
            omodel.params.(pn){1} = pv; pout{pp} = [pn ': ' pstr];
        end
        omodel = update_param(omodel, param_names{pp}, 0);
        save(fpath_opt, 'params', '-v7.3'); fprintf([strjoin(pout, ', ') '\n']);
    end
    % copy optimized parameters from session to group model objects
    for pp = 1:length(param_names)
        pn = param_names{pp}; model.params.(pn){ss} = omodel.params.(pn){1};
    end
end
model = update_param(model, param_names{pp}, 0);
model = pred_runs(model); model = pred_trials(model);
[roi, model] = tch_fit(roi, model, 0);

end
