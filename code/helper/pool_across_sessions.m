function roi_pool = pool_across_sessions(roi)
% Pools tchROI results across multiple scan sessions in each subject.
% 
% AS 5/2017

all_sessions = roi.sessions; subject_ids = {};
for ss = 1:length(all_sessions)
    [~, id] = fileparts(all_sessions{ss});
    subject_ids{ss} = id(1:3);
end

subjects = unique(subject_ids); subject_idxs = {};
for ss = 1:length(subjects)
    subject_idxs{ss} = find(strcmp(subjects{ss}, subject_ids));
end

sessions_init = cellfun(@(X) X(1), subject_idxs); vox_cnts = {}; vox_props = {};
roi_pool = tchROI(roi.name, roi.experiments, all_sessions(sessions_init));
for subject = 1:length(subjects)
    for session = 1:length(subject_idxs{subject})
        vox_cnts{subject}(session) = size(roi.runs{1, subject_idxs{subject}(session)}, 2);
    end
    vox_props{subject} = vox_cnts{subject} / sum(vox_cnts{subject});
end

% pool trial predictions and responses
for subject = 1:length(subjects)
    roi_pool.sessions{subject} = subjects{subject};
    for ee = 1:length(roi.experiments)
        empty_preds = cellfun(@(X) zeros(size(X)), roi.pred(:, subject_idxs{subject}(1), ee), 'uni', false);
        [predS, predT, predD, pred] = deal(empty_preds);
        trials = cellfun(@(X) zeros(size(X)), roi.trials(:, subject_idxs{subject}(1), ee), 'uni', false);
        for session = 1:length(subject_idxs{subject})
            if ~isempty(roi.predS)
                predS = cellfun(@(X, Y) X + Y * vox_props{subject}(session), predS, roi.predS(:, subject_idxs{subject}(session), ee), 'uni', false);
            end
            if ~isempty(roi.predT)
                predT = cellfun(@(X, Y) X + Y * vox_props{subject}(session), predT, roi.predT(:, subject_idxs{subject}(session), ee), 'uni', false);
            end
            if ~isempty(roi.predD)
                predD = cellfun(@(X, Y) X + Y * vox_props{subject}(session), predD, roi.predD(:, subject_idxs{subject}(session), ee), 'uni', false);
            end
            if ~isempty(roi.pred)
                pred = cellfun(@(X, Y) X + Y * vox_props{subject}(session), pred, roi.pred(:, subject_idxs{subject}(session), ee), 'uni', false);
            end
            if ~isempty(roi.trials)
                trials = cellfun(@(X, Y) X + Y * vox_props{subject}(session), trials, roi.trials(:, subject_idxs{subject}(session), ee), 'uni', false);
            end
        end
        if ~isempty(roi.predS)
            roi_pool.predS(:, subject, ee) = predS;
        end
        if ~isempty(roi.predT)
            roi_pool.predT(:, subject, ee) = predT;
        end
        if ~isempty(roi.predD)
            roi_pool.predD(:, subject, ee) = predD;
        end
        if ~isempty(roi.pred)
            roi_pool.pred(:, subject, ee) = pred;
        end
        if ~isempty(roi.trials)
            roi_pool.trials(:, subject, ee) = trials;
        end
    end
end

% pool model solutions and parameters
roi_pool.model.betas = {};
roi_pool.model.stdevs = {};
roi_pool.model.rbetas = {};
roi_pool.model.rstdevs = {};
roi_pool.model.varexp = {};
roi_pool.model.params = [];

roi_pool.model.type = roi.model.type;
roi_pool.model.cond_list = roi.model.cond_list;
roi_pool.model.cat_list = roi.model.cat_list;
roi_pool.model.pre_dur = roi.model.pre_dur;
roi_pool.model.post_dur = roi.model.post_dur;
roi_pool.model.fit_exps = roi.model.fit_exps;
roi_pool.model.num_channels = roi.model.num_channels;

for subject = 1:length(subjects)
    betas = zeros(size(roi.model.betas{subject_idxs{subject}(1)}));
    stdevs = zeros(size(roi.model.stdevs{subject_idxs{subject}(1)}));
    rbetas = zeros(size(roi.model.rbetas{subject_idxs{subject}(1)}));
    rstdevs = zeros(size(roi.model.rstdevs{subject_idxs{subject}(1)}));
    varexp = 0;
    param_names = fieldnames(roi.model.params);
    params = zeros(1, length(param_names));
    for session = 1:length(subject_idxs{subject})
        betas = betas + vox_props{subject}(session) * roi.model.betas{subject_idxs{subject}(session)};
        stdevs = stdevs + vox_props{subject}(session) * roi.model.stdevs{subject_idxs{subject}(session)};
        rbetas = rbetas + vox_props{subject}(session) * roi.model.rbetas{subject_idxs{subject}(session)};
        rstdevs = rstdevs + vox_props{subject}(session) * roi.model.rstdevs{subject_idxs{subject}(session)};
        varexp = varexp + vox_props{subject}(session) * roi.model.varexp{subject_idxs{subject}(session)};
        for pp = 1:length(param_names)
            params(pp) = params(pp) + vox_props{subject}(session) * roi.model.params.(param_names{pp}){subject_idxs{subject}(session)};
        end
    end
    roi_pool.model.betas{subject} = betas;
    roi_pool.model.stdevs{subject} = stdevs;
    roi_pool.model.rbetas{subject} = rbetas;
    roi_pool.model.rstdevs{subject} = rstdevs;
    roi_pool.model.varexp{subject} = varexp;
    roi_pool.model.params = params;
end

end
