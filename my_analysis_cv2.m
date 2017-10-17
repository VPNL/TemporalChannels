function my_analysis_cv2(rois, models, sessions)

if nargin < 1 || isempty(rois)
    rois = {'IOG_faces' 'pFus_faces' 'mFus_faces' 'pSTS_faces'...
        'LOS_bodies' 'ITG_bodies' 'MTG_bodies' 'OTS_bodies' ...
        'IOS_characters' 'pOTS_characters' 'mOTS_characters' ...
        'V1' 'V2' 'V3' 'hV4' 'MT'};
else
    rois = force_cell(rois);
end

if nargin < 2 || isempty(models)
    models = {'1ch-lin' '1ch-pow' '1ch-exp' ...
        '2ch-lin-quad' '2ch-pow-quad' '2ch-exp-quad' ...
        '3ch-lin-quad-exp' '3ch-pow-quad-exp' '3ch-exp-quad-exp'};
else
    models = force_cell(models);
end


if nargin < 3
    sessions = {'as' 'bj' 'cs' 'md' 'em' 'jg' 'kg' 'md' 'mg' 'wd' 'yl'};
else
    force_cell(sessions);
end

fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
val_exps = {'ExpAe' 'ExpBe' 'ExpCe'; 'ExpAo' 'ExpBo' 'ExpCo'};

project_dir = fullfile(RAID, 'projects', 'CategoryChannels', 'model');
for rr = 1:length(rois)
    for mm = 1:length(models)
        roi = tchROI(rois{rr}, fit_exps(1, :));
        fname = [roi(1).nickname '_' models{mm} '_split-half-concat.mat'];
        fpath = fullfile(project_dir, 'TemporalChannels', 'results', fname);
        if exist(fpath, 'file') == 0
            [roi1, model1] = tch_model_roi(rois{rr}, models{mm}, ...
                fit_exps(1, :), val_exps(1, :), sessions);
            [roi2, model2] = tch_model_roi(rois{rr}, models{mm}, ...
                fit_exps(2, :), val_exps(2, :), sessions);
            roi = pool_across_folds(roi1, roi2);
            save(fpath, 'roi', '-v7.3');
        end
    end
end

end
