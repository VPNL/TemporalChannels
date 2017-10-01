function make_model_grids(model_type, exps)

roi = ROI('pFus_faces', exps);
roi = select_sessions(roi);
sessions = roi.sessions;

for ss = 1:length(sessions)
    roi = ROI('pFus_faces', exps, sessions{ss});
    roi = tc_runs(roi);
    model = ModelTS(model_type, roi.experiments, sessions{ss});
    model = code_stim(model);
    model = pred_runs(model);
    model = pred_trials(model);
    roi = tc_trials(roi, model);
    [roi, model] = tc_fit(roi, model);
    [rois, models] = grid_search(roi, model, 1, 5);
end

end
