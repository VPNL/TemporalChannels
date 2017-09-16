function model = pred_trials_balloon(model)
% Generates trial predictors using the standard model.

% define global variables and signal parameters
global which_tau tau_n tau_p tau tauMTT alpha E0 V0;
% get design parameters
sessions = model.sessions; nsess = length(sessions);
params = model.params; irfs = model.irfs; dt = params.delta_t;
fs = model.fs; tr = model.tr; cond_list = model.cond_list;
stimfiles = model.stimfiles; nruns = model.num_runs;
model.trial_preds.pred = cell(max(cellfun(@length, cond_list)), nsess, model.num_exps);

% generate trial predictors for each session
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
        t = 0:params.delta_t:model.pre_dur + td + model.post_dur;
        for ss = 1%:length(sessions)
            in_flow = convolve_vecs(cstim, irfs.gamma{ss}, fs, fs);
            in_flow = 0.7 * (in_flow / sum(irfs.gamma{ss})) + 1;
            for pp = 1:size(cstim, 2)
                % initialize variables
                [v, q, IN_FLOW, OUT_FLOW, CMRO2] = deal(1);
                S = 0; OEF = params.E0;
                which_tau = 1; tau = params.tau_p; tauMTT = params.tauMTT;
                tau_p = params.tau_p; tau_n = params.tau_n;
                E0 = params.E0; V0 = params.V0; alpha = params.alpha;
                % get the simulated values of all variables                
                for ii = 1:length(t) - 1
                    ii_flow = in_flow(ii, pp);
                    v(ii + 1) = runge_kutta(dt, @dvdt, t(ii), v(ii), ii_flow);
                    OUT_FLOW(ii + 1) = flow_out(v(ii + 1), t(ii), ii_flow);
                    q(ii + 1) = runge_kutta(dt, @dqdt, t(ii), q(ii), v(ii), ii_flow);  
                    S1 = params.k1 * (1 - q(ii + 1));
                    S2 = params.k2 * (1 - q(ii + 1) / v(ii + 1));
                    S3 = params.k3 * (1 - v(ii + 1));
                    S(ii + 1) = V0 * (S1 + S2 + S3);
                    OEF(ii + 1) = 1 - (1 - E0) .^ (1 ./ ii_flow);
                    CMRO2(ii + 1) = (OEF(ii + 1) / E0) * IN_FLOW(ii);
                    IN_FLOW(ii + 1) = ii_flow;
                end
                Sr = convolve_vecs(S(1:length(t) - 1)', 1, fs, 1 / tr);
                model.trial_preds.pred{cc, ss, ee} = Sr;
            end
        end
        for ss = 2:nsess
            model.trial_preds.pred{cc, ss, ee} = model.trial_preds.pred{cc, 1, ee};
        end
    end
    rcnt = rcnt + nruns(ee, 1);
end

end
