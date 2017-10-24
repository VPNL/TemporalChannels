[mdir, ~] = fileparts(pwd); addpath(genpath(mdir));
model_type = '2ch-pow-rect';
rois = {'IOG_faces' 'pFus_faces' 'mFus_faces' 'pSTS_faces' ...
    'LOS_bodies' 'ITG_bodies' 'MTG_bodies' 'OTS_bodies' ...
    'IOS_characters' 'pOTS_characters' 'mOTS_characters'};
fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
val_exps = {'ExpAe' 'ExpBe' 'ExpCe'; 'ExpAo' 'ExpBo' 'ExpCo'};
sessions = {'as' 'bj' 'cs' 'em' 'jg' 'kg' 'md' 'mg' 'sc' 'wd' 'yl'};

for rr = 1:length(rois)
    fname = [strrep(rois{rr}, '_', '-') '_' model_type '_split-half.mat'];
    if exist(fullfile(mdir, 'results', fname), 'file') == 0
        [roi1, model1] = tch_model_roi(rois{rr}, model_type, ...
            fit_exps(1, :), val_exps(1, :)', 1, sessions);
        [roi2, model2] = tch_model_roi(rois{rr}, model_type, ...
            fit_exps(2, :), val_exps(2, :)', 1, sessions);
        roi = pool_across_folds(roi1, roi2);
        save(fullfile(mdir, 'results', fname), 'roi');
    end
end
