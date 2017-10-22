function [roi, model] = tch_optimize_grid_grad(roi, model, fit_exps)
% Custom two-stage optimization prodcedure consisting of grid search
% followed by iterative gradident descent from multiple seed points. 

param_names = fieldnames(model.params); sessions = roi.sessions;
grid_stem = 'grid_search_results_'; grad_stem = 'grad_desc_results_';
for ss = 1:length(sessions)
    fname_grid = [grid_stem model.type '_fit' [fit_exps{:}] '.mat'];
    fpath_grid = fullfile(sessions{ss}, 'ROIs', roi.name, fname_grid);
    fname_grad = [grad_stem model.type '_fit' [fit_exps{:}] '.mat'];
    fpath_grad = fullfile(sessions{ss}, 'ROIs', roi.name, fname_grad);
    % load optimization results if saved, otherwise compute
    if exist(fpath_grad, 'file') == 2
        load(fpath_grad); fprintf('Loading gradient descent results. \n');
        sroi = tch_runs(tchROI(roi.name, roi.experiments, sessions{ss}));
        smodel = tchModel(model.type, roi.experiments, sessions{ss});
        smodel.normT = model.normT; smodel.normD = model.normD;
        for pp = 1:length(param_names)
            pn = param_names{pp}; smodel.params.(pn){1} = params.(pn){1};
        end
        smodel = code_stim(smodel); smodel = update_param(smodel, pn, 0);
        smodel = pred_runs(smodel); smodel = pred_trials(smodel);
    elseif exist(fpath_grid, 'file') == 2
        load(fpath_grid); fprintf('Loading grid search results. \n');
        sroi = tch_runs(tchROI(roi.name, roi.experiments, sessions{ss}));
        smodel = tchModel(model.type, roi.experiments, sessions{ss});
        smodel.normT = model.normT; smodel.normD = model.normD;
        smodel = code_stim(smodel);
        for mm = 1:length(params)
            smodel(mm) = smodel(1); sroi(mm) = sroi(1);
            for pp = 1:length(param_names)
                pn = param_names{pp};
                smodel(mm).params.(pn){1} = params(mm).(pn){1};
            end
            smodel(mm) = update_param(smodel(mm), pn, 0);
            smodel(mm) = pred_runs(smodel(mm));
            smodel(mm) = pred_trials(smodel(mm));
            sroi(mm) = tch_trials(sroi(mm), smodel(mm));
            sroi(mm) = tch_fit(sroi(mm), smodel(mm));
        end
        [~, smodel] = tch_optimize_fit(sroi, smodel);
        params = smodel(1).params; save(fpath_grad, 'params', '-v7.3');
    else
        [sroi, smodel] = tch_grid_search(roi, model, ss, 5);
        params = smodel(1).params;
        for mm = 2:length(smodel); params(mm) = smodel(mm).params; end
        save(fpath_grid, 'params', '-v7.3');
        [~, smodel] = tch_optimize_fit(sroi, smodel);
        params = smodel(1).params; save(fpath_grad, 'params', '-v7.3');
    end
    % copy optimized parameters from session to group model objects
    for pp = 1:length(param_names)
        model.params.(param_names{pp}){ss} = smodel.params.(param_names{pp}){1};
    end
end
model = update_param(model, param_names{pp}, 0); model = norm_model(model, 1);
model = pred_runs(model); model = pred_trials(model);
[roi, model] = tch_fit(roi, model, 0);

end
