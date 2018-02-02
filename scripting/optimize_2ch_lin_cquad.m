[mdir, ~] = fileparts(pwd); addpath(genpath(mdir));
model_type = '2ch-lin-cquad';
rois = {'mFus_faces' 'OTS_bodies' 'mOTS_characters' ...
    'pSTS_faces' 'MTG_bodies' 'V1' 'hV4' 'MT' ...
    'pFus_faces' 'ITG_bodies' 'pOTS_characters' ...
    'IOG_faces' 'LOS_bodies' 'IOS_characters'};
fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
val_exps = {'ExpAe' 'ExpBe' 'ExpCe'; 'ExpAo' 'ExpBo' 'ExpCo'};
sessions = {'ad' 'as' 'bj' 'cs' 'em' 'jg' 'kg' 'md' 'mg' 'sc' 'wd' 'yl'};

for rr = 1:length(rois)
    fname = [strrep(rois{rr}, '_', '-') '_' model_type '_split-half.mat'];
    if exist(fullfile(mdir, 'results', fname), 'file') == 0
        [roi1, model1] = tch_model_roi(rois{rr}, model_type, ...
            fit_exps(1, :), val_exps(1, :), 1, sessions);
        [roi2, model2] = tch_model_roi(rois{rr}, model_type, ...
            fit_exps(2, :), val_exps(2, :), 1, sessions);
        roi = pool_across_folds(roi1, roi2);
        save(fullfile(mdir, 'results', fname), 'roi');
    end
end
