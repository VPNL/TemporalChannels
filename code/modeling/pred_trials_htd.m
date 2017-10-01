function model = pred_trials_htd(model)
% Generates trial predictors using the hemodynamic temporal derivative 
% (HTD) model proposed by Henson et al. (2002). 

% get design parameters
sessions = model.sessions; nsess = length(sessions);
params = model.params; irfs = model.irfs;
fs = model.fs; tr = model.tr; cond_list = model.cond_list;
stimfiles = model.stimfiles; nruns = model.num_runs;
model.trial_preds.S = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);
model.trial_preds.T = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);

rcnt = 1;
for ee = 1:model.num_exps
    [on, off, c, ims, ton, toff, tc, rd, cl] = tch_stimfile(stimfiles{rcnt, 1});
    istim = model.stim{rcnt, 1};
    for cc = 1:length(cond_list{ee})
        % find trial onset and offset times
        ii = find(strcmp(cl{cc}, tc), 1);
        ion = ton(ii); ioff = ceil(toff(ii) - .001); td = ioff - ion;
        % extract stimulus vector from condition time window
        cstim = istim(fs * (ion - model.pre_dur) + 1:round(fs * (ion + td + model.post_dur)), :);
        for ss = 1:length(sessions)
            % convolve stimulus with hrf and dhrf
            fmri_c = convolve_vecs(cstim,  irfs.hrf{ss}, fs, 1 / tr);
            fmri_d = convolve_vecs(cstim,  irfs.dhrf{ss}, fs, 1 / tr);
            % store fMRI predictors in model structure
            model.trial_preds.S{cc, ss, ee} = fmri_c;
            model.trial_preds.T{cc, ss, ee} = fmri_d;
        end
    end
    rcnt = rcnt + nruns(ee, 1);
end

end
