function model = trial_preds_latency(model)
% Generates trial predictors using the 2 temporal-channel model.

% get design parameters
sessions = model.sessions; nsess = length(sessions);
params = model.params; irfs = model.irfs;
fs = model.fs; tr = model.tr; cond_list = model.cond_list;
stimfiles = model.stimfiles; nruns = model.num_runs;
model.tc_pred.S = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);
model.tc_pred.T = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);

rcnt = 1;
for ee = 1:model.num_exps
    [on, off, c, ims, ton, toff, tc, rd, cl] = stimfileTS(stimfiles{rcnt, 1});
    istim = model.stim{rcnt, 1};
    for cc = 1:length(cond_list{ee})
        % find trial onset and offset times
        ii = find(strcmp(cl{cc}, tc), 1);
        ion = ton(ii); ioff = ceil(toff(ii) - .001); td = ioff - ion;
        % extract stimulus vector from condition time window
        cstim = istim(fs * (ion - model.pre_dur) + 1:round(fs * (ion + td + model.post_dur)), :);
        for ss = 1:length(sessions)
            % convolve stimulus with hrf and dhrf
            fmri_c = convolve_vecs(cstim,  irfs.hrf{ss}, fs, 1 / model.tr);
            fmri_d = convolve_vecs(cstim,  irfs.dhrf{ss}, fs, 1 / model.tr);
            % store fMRI predictors in model structure
            model.tc_pred.S{cc, ss, ee} = fmri_c;
            model.tc_pred.T{cc, ss, ee} = fmri_d;
        end
    end
    rcnt = rcnt + nruns(ee, 1);
end

end
