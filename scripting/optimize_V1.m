[mdir, ~] = fileparts(pwd); addpath(genpath(mdir));
models = {'3ch-exp-rect-exp' '3ch-pow-rect-exp' '3ch-lin-rect-exp' ...
    '2ch-exp-rect' '2ch-pow-rect' '1ch-pow' '1ch-exp'};
rois = {'V1' 'hV4' 'MT'};
fit_exps = {'ExpAo' 'ExpBo' 'ExpCo'; 'ExpAe' 'ExpBe' 'ExpCe'};
val_exps = {'ExpAe' 'ExpBe' 'ExpCe'; 'ExpAo' 'ExpBo' 'ExpCo'};
sessions = {'as' 'bj' 'cs' 'em' 'jg' 'kg' 'md' 'mg' 'sc' 'wd' 'yl'};

for rr = 1:length(rois)
    for mm = 1:length(models)
        fname = [strrep(rois{rr}, '_', '-') '_' models{mm} '_split-half.mat'];
        if exist(fullfile(mdir, 'results', fname), 'file') == 0
            [roi1, model1] = tch_model_roi(rois{rr}, models{mm}, ...
                fit_exps(1, :), val_exps(1, :)', 1, sessions);
            [roi2, model2] = tch_model_roi(rois{rr}, models{mm}, ...
                fit_exps(2, :), val_exps(2, :)', 1, sessions);
            roi = pool_across_folds(roi1, roi2);
            save(fullfile(mdir, 'results', fname), 'roi');
        end
    end
end
