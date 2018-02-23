function [roi, model] = tch_optimize_fmincon(roi, model, fit_exps)
% Nonlinear optimization prodcedure using built-in fmincon algorithm.

param_names = fieldnames(model.params); sessions = roi.sessions;
if length(param_names) <= 3; fprec = 1e-8; end
if length(param_names) > 3; fprec = 1e-6; end
for ss = 1:length(sessions)
    fname_opt = ['optimization_results_' model.type '_fit' [fit_exps{:}] '.mat'];
    fpath_opt = fullfile(sessions{ss}, 'ROIs', roi.name, fname_opt);
    % load optimization results if saved, otherwise compute
    if exist(fpath_opt, 'file') == 2
        load(fpath_opt); fprintf('Loading gradient descent results. \n');
        omodel = tchModel(model.type, roi.experiments, sessions{ss});
        omodel.tr = model.tr; omodel = code_stim(omodel);
        for pp = 1:length(param_names)
            pn = param_names{pp}; omodel.params.(pn){1} = params.(pn){1};
        end
        omodel = tch_update_param(omodel, pn, 0);
    else
        sroi = tchROI(roi.name, roi.experiments, sessions{ss});
        sroi.tr = roi.tr; sroi = tch_runs(sroi);
        omodel = tchModel(model.type, roi.experiments, sessions{ss});
        omodel.tr = model.tr; omodel = code_stim(omodel);
        omodel.normT = model.normT; omodel.normP = model.normP;
        omodel = pred_runs(omodel); omodel = pred_trials(omodel);
        sroi = tch_trials(sroi, omodel); sroi = tch_fit(sroi, omodel);
        pout = cell(1, length(param_names));
        fmin_options = optimoptions('fmincon', 'Display', 'off', ...
            'StepTolerance', fprec, 'UseParallel', true);
        fprintf('Optimizing parameters for %s ...\n', roi.session_ids{ss});
        switch model.type
            case '1ch-lin'
                obj_fun = tch_obj_fun_1ch_lin(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                x_init = tau_s; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_s{1} = x_opt(1);
            case '1ch-exp'
                obj_fun = tch_obj_fun_1ch_exp(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_ae]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1], [20 4], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 10000;
            case '1ch-rect'
                obj_fun = tch_obj_fun_1ch_rect(sroi, omodel);
                x_init = omodel.params.tau_t{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_t{1} = x_opt(1);
            case '1ch-quad'
                obj_fun = tch_obj_fun_1ch_quad(sroi, omodel);
                x_init = omodel.params.tau_t{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_t{1} = x_opt(1);
            case '1ch-cquad'
                obj_fun = tch_obj_fun_1ch_cquad(sroi, omodel);
                x_init = omodel.params.tau_t{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_t{1} = x_opt(1);
            case '1ch-sig'
                obj_fun = tch_obj_fun_1ch_sig(sroi, omodel);
                tau_t = omodel.params.tau_t{1};
                lambda_p = omodel.params.lambda_p{1};
                kappa_p = omodel.params.kappa_p{1};
                kappa_n = omodel.params.kappa_n{1};
                x_init = [tau_t lambda_p kappa_p kappa_n];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 .01 .1 .1], [20 .5 6 6], [], fmin_options);
                params.tau_t{1} = x_opt(1);
                params.lambda_p{1} = x_opt(2);
                params.kappa_p{1} = x_opt(3);
                params.kappa_n{1} = x_opt(4);
            case '2ch-lin-rect'
                obj_fun = tch_obj_fun_2ch_lin_rect(sroi, omodel);
                x_init = omodel.params.tau_s{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_s{1} = x_opt(1);
            case '2ch-exp-rect'
                obj_fun = tch_obj_fun_2ch_exp_rect(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_ae]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1], [20 4], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 10000;
            case '2ch-lin-quad'
                obj_fun = tch_obj_fun_2ch_lin_quad(sroi, omodel);
                x_init = omodel.params.tau_s{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_s{1} = x_opt(1);
            case '2ch-exp-quad'
                obj_fun = tch_obj_fun_2ch_exp_quad(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_ae]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1], [20 4], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 10000;
            case '2ch-lin-cquad'
                obj_fun = tch_obj_fun_2ch_lin_cquad(sroi, omodel);
                x_init = omodel.params.tau_s{1};
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], 4, 20, [], fmin_options);
                params.tau_s{1} = x_opt(1);
            case '2ch-exp-cquad'
                obj_fun = tch_obj_fun_2ch_exp_cquad(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_ae];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1], [20 4], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 10000;
            case '2ch-lin-sig'
                obj_fun = tch_obj_fun_2ch_lin_sig(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                lambda_p = omodel.params.lambda_p{1};
                kappa_p = omodel.params.kappa_p{1};
                kappa_n = omodel.params.kappa_n{1};
                x_init = [tau_s lambda_p kappa_p kappa_n];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 .01 .1 .1], [20 .5 6 6], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.lambda_p{1} = x_opt(2);
                params.kappa_p{1} = x_opt(3);
                params.kappa_n{1} = x_opt(4);
            case '2ch-exp-sig'
                obj_fun = tch_obj_fun_2ch_exp_sig(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                lambda_p = omodel.params.lambda_p{1};
                kappa_p = omodel.params.kappa_p{1};
                kappa_n = omodel.params.kappa_n{1};
                x_init = [tau_s tau_ae lambda_p kappa_p kappa_n];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1 .01 .1 .1], [20 4 .5 6 6], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 10000;
                params.lambda_p{1} = x_opt(3);
                params.kappa_p{1} = x_opt(4);
                params.kappa_n{1} = x_opt(5);
            case '2ch-exp-crect'
                obj_fun = tch_obj_fun_2ch_exp_crect(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_t = omodel.params.tau_t{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_t tau_ae];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1 1], [20 20 3], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_t{1} = x_opt(2);
                params.tau_ae{1} = x_opt(3) * 10000;
            case '3ch-exp-crect-crect'
                obj_fun = tch_obj_fun_3ch_exp_crect_crect(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_t = omodel.params.tau_t{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_t tau_ae];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1 1], [20 20 4], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_t{1} = x_opt(2);
                params.tau_ae{1} = x_opt(3) * 10000;
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
            case '3ch-lin-quad-exp'
                obj_fun = tch_obj_fun_3ch_lin_quad_exp(sroi, omodel);
                x_init = omodel.params.tau_pe{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], .1, 12, [], fmin_options);
                params.tau_pe{1} = x_opt(1) * 1000;
            case '3ch-lin-rect-exp'
                obj_fun = tch_obj_fun_3ch_lin_rect_exp(sroi, omodel);
                x_init = omodel.params.tau_pe{1} / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], .1, 12, [], fmin_options);
                params.tau_pe{1} = x_opt(1) * 1000;
            case '3ch-pow-quad-exp'
                obj_fun = tch_obj_fun_3ch_pow_quad_exp(sroi, omodel);
                x_init = [omodel.params.epsilon{1} omodel.params.tau_pe{1} / 1000];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [.001 .1], [1 12], [], fmin_options);
                params.epsilon{1} = x_opt(1);
                params.tau_pe{1} = x_opt(2) * 1000;
            case '3ch-pow-rect-exp'
                obj_fun = tch_obj_fun_3ch_pow_rect_exp(sroi, omodel);
                x_init = [omodel.params.epsilon{1} omodel.params.tau_pe{1} / 1000];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [.001 .1], [1 12], [], fmin_options);
                params.epsilon{1} = x_opt(1);
                params.tau_pe{1} = x_opt(2) * 1000;
             case '3ch-exp-quad-exp'
                obj_fun = tch_obj_fun_3ch_exp_quad_exp(sroi, omodel);
                x_init = [omodel.params.tau_ae{1} omodel.params.tau_pe{1}] / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [1 .1], [60 12], [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
                params.tau_pe{1} = x_opt(2) * 1000;
            case '3ch-exp-rect-exp'
                obj_fun = tch_obj_fun_3ch_exp_rect_exp(sroi, omodel);
                x_init = [omodel.params.tau_ae{1} omodel.params.tau_pe{1}] / 1000;
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [1 .1], [60 12], [], fmin_options);
                params.tau_ae{1} = x_opt(1) * 1000;
                params.tau_pe{1} = x_opt(2) * 1000;
            case '3ch-lin-quad-exp-opt'
                obj_fun = tch_obj_fun_3ch_lin_quad_exp_opt(sroi, omodel);
                tau_s = omodel.params.tau_s{1}; kappa = omodel.params.kappa{1};
                n1 = omodel.params.n1{1}; n2 = omodel.params.n2{1};
                tau_pe = omodel.params.tau_pe{1} / 1000;
                x_init = [tau_s n1 n2 kappa tau_pe]; 
                x_opt = fmincon(obj_fun, x_init, [0 1 -1 0 0], 0, ...
                    [], [], [4 9 10 1 .1], [50 19 20 10 12], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.n1{1} = x_opt(2);
                params.n2{1} = x_opt(3);
                params.kappa{1} = x_opt(4);
                params.tau_pe{1} = x_opt(5) * 1000;
            case '3ch-lin-rect-exp-opt'
                obj_fun = tch_obj_fun_3ch_lin_rect_exp_opt(sroi, omodel);
                tau_s = omodel.params.tau_s{1}; kappa = omodel.params.kappa{1};
                n1 = omodel.params.n1{1}; n2 = omodel.params.n2{1};
                tau_pe = omodel.params.tau_pe{1} / 1000;
                x_init = [tau_s n1 n2 kappa tau_pe]; 
                x_opt = fmincon(obj_fun, x_init, [0 1 -1 0 0], 0, ...
                    [], [], [4 9 10 1 .1], [50 19 20 10 12], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.n1{1} = x_opt(2);
                params.n2{1} = x_opt(3);
                params.kappa{1} = x_opt(4);
                params.tau_pe{1} = x_opt(5) * 1000;
            case '3ch-exp-quad-exp-opt'
                obj_fun = tch_obj_fun_3ch_exp_quad_exp_opt(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 1000;
                tau_pe = omodel.params.tau_pe{1} / 1000;
                x_init = [tau_s tau_ae tau_pe]; 
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1 .1], [20 60 12], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 1000;
                params.tau_pe{1} = x_opt(3) * 1000;
            case '3ch-exp-quad-crect-opt'
                obj_fun = tch_obj_fun_3ch_exp_quad_crect_opt(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 1000;
                tau_p = omodel.params.tau_p{1};
                epsilon = omodel.params.epsilon{1};
                x_init = [tau_s tau_ae tau_p epsilon];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 12 4 .001], [20 60 20 1], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 1000;
                params.tau_p{1} = x_opt(3);
                params.epsilon{1} = x_opt(4);
            case '2ch-lin-crect'
                obj_fun = tch_obj_fun_2ch_lin_crect(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                epsilon = omodel.params.epsilon{1};
                x_init = [tau_s epsilon];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 .001], [20 .5], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.epsilon{1} = x_opt(2);
            case '2ch-exp-dquad'
                obj_fun = tch_obj_fun_2ch_exp_dquad(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_ae = omodel.params.tau_ae{1} / 1000;
                sigma = omodel.params.sigma{1};
                x_init = [tau_s tau_ae sigma];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 12 .001], [20 60 .5], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_ae{1} = x_opt(2) * 1000;
                params.sigma{1} = x_opt(3);
            case '2ch-lin-dquad'
                obj_fun = tch_obj_fun_2ch_lin_dquad(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                sigma = omodel.params.sigma{1};
                x_init = [tau_s sigma];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 .001], [20 .5], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.sigma{1} = x_opt(2);
            case '3ch-exp-cquad-rect'
                obj_fun = tch_obj_fun_3ch_exp_cquad_rect(sroi, omodel);
                tau_s = omodel.params.tau_s{1};
                tau_t = omodel.params.tau_t{1};
                tau_p = omodel.params.tau_p{1};
                tau_ae = omodel.params.tau_ae{1} / 10000;
                x_init = [tau_s tau_t tau_p tau_ae];
                x_opt = fmincon(obj_fun, x_init, [], [], ...
                    [], [], [4 1 8 1], [20 20 24 6], [], fmin_options);
                params.tau_s{1} = x_opt(1);
                params.tau_t{1} = x_opt(2);
                params.tau_p{1} = x_opt(3);
                params.tau_ae{1} = x_opt(4) * 10000;
        end
        for pp = 1:length(param_names)
            pn = param_names{pp}; pv = params.(pn){1}; pstr = num2str(pv, 3);
            omodel.params.(pn){1} = pv; pout{pp} = [pn ': ' pstr];
        end
        omodel = tch_update_param(omodel, param_names{pp}, 0);
        save(fpath_opt, 'params', '-v7.3'); fprintf([strjoin(pout, ', ') '\n']);
    end
    % copy optimized parameters from session to group model objects
    for pp = 1:length(param_names)
        pn = param_names{pp}; model.params.(pn){ss} = omodel.params.(pn){1};
    end
end
model = tch_update_param(model, param_names{pp}, 0);
model = pred_runs(model); model = pred_trials(model);
[roi, model] = tch_fit(roi, model, 0);

end
