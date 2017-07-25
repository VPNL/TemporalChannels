function model = trial_preds_all(model)
% Generates trial predictors using the ______ model. 

% get design parameters
sessions = model.sessions; nsess = length(sessions);
params = model.params; irfs = model.irfs;
fs = model.fs; tr = model.tr; cond_list = model.cond_list;
stimfiles = model.stimfiles; nruns = model.num_runs;

if sum(strcmp(model.type, {'standard' 'cts' 'dcts'}))
    model.tc_pred.pred = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);
else
    model.tc_pred.S = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);
    model.tc_pred.T = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);
    if strcmp(model.type, '3ch')
        model.tc_pred.D = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);
    end
end

rcnt = 1;
for ee = 1:model.num_exps
    [on, off, c, ims, ton, toff, tc, rd, cl] = stimfileTS(stimfiles{rcnt, 1});
    istim = model.stim{rcnt, 1};
    for cc = 1:length(cond_list{ee})
        % find trial onset and offset times
        ii = find(strcmp(cl{cc}, tc), 1);
        ion = ton(ii); ioff = ceil(toff(ii)-.001); td = ioff-ion;
        % extract stimulus vector from condition time window
        cstim = istim(fs * (ion - model.pre_dur) + 1:round(fs * (ion + td + model.post_dur)), :);
        for ss = 1:length(sessions)
            % if using 2ch model
            if strcmp(model.type, 'standard')
                % if using standard model convolve stimulus with HRF
                fmri = convolve_vecs(cstim, fs, model.tr, irfs.hrf{ss});
                % store fMRI predictor in model structure
                model.tc_pred.pred{cc, ss, ee} = fmri;
            elseif sum(strcmp(model.type, {'2ch' '2ch-cts' '2ch-dcts' '3ch'}))
                % convolve stimulus with channel IRFs
                predS = convolve_vecs(cstim, fs, 1/fs, irfs.nrfS{ss});
                predT = convolve_vecs(cstim, fs, 1/fs, irfs.nrfT{ss}) .^ 2;
                switch model.type
                    case '2ch-cts'
                        predS = predS.^params.e{ss};
                    case '2ch-dcts'
                        predSn = predS.^2;
                        predSd = convolve_vecs(predS, fs, 1/fs, irfs.lpf{ss});
                        predSd = params.sigma{ss}.^2 + predSd.^2;
                        predS = predSn./predSd;
                    case '3ch'
                        dstim = double(~cstim);
                        predD = convolve_vecs(dstim, fs, 1/fs, irfs.nrfD{ss});
                        predDn = predD.^2;
                        predDd = convolve_vecs(predD, fs, 1/fs, irfs.lpf{ss});
                        predDd = params.sigma{ss}.^2 + predDd.^2;
                        predD = predDn./predDd;
                        fmriD = convolve_vecs(predD, fs, tr, irfs.hrf{ss});
                        model.tc_pred.D{cc, ss, ee} = fmriD*model.normD;
                        predS = predS.^params.e{ss};
                end
                % convolve neural predictors with HRF
                fmriS = convolve_vecs(predS, fs, tr, irfs.hrf{ss});
                fmriT = convolve_vecs(predT, fs, tr, irfs.hrf{ss});
                % store fMRI predictors in model structure
                model.tc_pred.S{cc, ss, ee} = fmriS;
                model.tc_pred.T{cc, ss, ee} = fmriT * model.normT;
            elseif strcmp(model.type, 'cts')
                % convolve stimulus with channel IRF and HRF
                predS = convolve_vecs(cstim, fs, 1/fs, irfs.nrfS{ss}).^params.e{ss};
                fmriS = convolve_vecs(predS, fs, tr, irfs.hrf{ss});
                model.tc_pred.pred{cc, ss, ee} = fmriS;
            elseif strcmp(model.type, 'dcts')
                % convolve stimulus with channel IRF and lpf
                predS = convolve_vecs(cstim, fs, 1/fs, irfs.nrfS{ss});
                predSn = predS.^2;
                predSd = convolve_vecs(predS, fs, 1/fs, irfs.lpf{ss});
                predSd = params.sigma{ss}.^2 + predSd.^2;
                predS = predSn./predSd;
                % convolve neural predictors with HRF
                fmriS = convolve_vecs(predS, fs, tr, irfs.hrf{ss});
                model.tc_pred.pred{cc, ss, ee} = fmriS;
            end
        end
    end
    rcnt = rcnt+nruns(ee, 1);
end

end
