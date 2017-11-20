roi_names = {'IOG_faces' 'pFus_faces' 'mFus_faces' ...
    'LOS_bodies' 'ITG_bodies' 'MTG_bodies' 'OTS_bodies' ...
    'IOS_characters' 'IOS_characters' 'IOS_characters'};
model_types = {'1ch-lin' '2ch-lin-quad' '2ch-lin-rect'};
exps = {'FastLocs_1Hz' 'FastLocs_2Hz' 'FastLocs_4Hz' 'FastLocs_8Hz'};

for rr = 1:length(roi_names)
    for mm = 1:length(model_types)
        roi = tchROI(roi_names{rr}, exps); roi.tr = 2;
        roi = tch_runs(roi);
        model = tchModel(model_types{mm}, roi.experiments, roi.sessions);
        model.tr = 2;
        model = code_stim(model);
        %model = norm_model(model, 1);
        model = pred_runs(model);
        model = pred_trials(model);
        roi = tch_trials(roi, model);
        roi = tch_fit(roi, model);
        roi = tch_pred(roi, model);
        fname = [roi(1).nickname '_' roi(1).model.type '_fit' [roi(1).experiments{:}]];
        save(fullfile(roi.project_dir, 'results', fname), 'roi', 'model', '-v7.3');
    end
end
