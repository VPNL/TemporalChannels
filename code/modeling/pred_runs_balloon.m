function model = pred_runs_balloon(model)

% define global variables
global which_tau tau_n tau_p tau tauMTT alpha E0 V0;
% get design parameters
params = model.params; irfs_init = model.irfs; dt = params.delta_t;
fs = model.fs; tr = model.tr; rd = model.run_durs; stim = model.stim;
cat_list = unique([model.cats{:}]); ncats = length(cat_list);
[nruns_max, nsess] = size(model.run_durs); empty_cells = cellfun(@isempty, rd);
irfs_names = fieldnames(irfs_init); irfs = struct;
for ff = 1:length(irfs_names)
    irfs.(irfs_names{ff}) = repmat(irfs_init.(irfs_names{ff}), nruns_max, 1);
    irfs.(irfs_names{ff})(empty_cells) = {[]};
end
time_vecs = cellfun(@(X) [0:dt:X]', rd, 'uni', false);
time_vecs(empty_cells) = {[]};

%% model nonlinear hemodynamics

% convolve input with a gamma and normalize
in_flow = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), stim, irfs.gamma, 'uni', false);
in_flow = cellfun(@(X, Y) 0.7 * (X / sum(Y)) + 1, in_flow, irfs.gamma, 'uni', false);
run_preds = cellfun(@(X) zeros(X / tr, ncats), rd, 'uni', false);
run_preds(empty_cells) = {[]};

% get the simulated values of all variables for each predictor in each run
for ss = 1:nsess
    [~, session_id] = fileparts(model.sessions{ss});
    fprintf('Simulating balloon model for %s', session_id);
    fname = ['balloon_model_' [model.experiments{:}] '.mat'];
    fpath = fullfile(model.sessions{ss}, 'Stimuli', fname);
    if exist(fpath, 'file') == 2
        load(fpath);
        run_preds(:, ss) = session_run_preds;
    else
        for rr = 1:nruns_max
            fprintf('.');
            for pp = 1:size(in_flow{rr, ss}, 2)
                % initialize variables
                [v, q, IN_FLOW, OUT_FLOW, CMRO2] = deal(1);
                S = 0; OEF = params.E0; t = time_vecs{rr, ss}; 
                which_tau = 1; tau = params.tau_i; tauMTT = params.tauMTT;
                tau_p = params.tau_p; tau_n = params.tau_n;
                E0 = params.E0; V0 = params.V0; alpha = params.alpha;
                % get the simulated values of all variables
                for ii = 1:size(time_vecs{rr, ss}, 1) - 1
                    ii_flow = in_flow{rr, ss}(ii, pp);
                    v(ii + 1) = runge_kutta(dt, @dvdt, t(ii), v(ii), ii_flow);
                    OUT_FLOW(ii + 1) = flow_out(v(ii + 1), t(ii), ii_flow);
                    q(ii + 1) = runge_kutta(dt, @dqdt, t(ii), q(ii), v(ii), ii_flow);  
                    S1 = params.k1 * (1 - q(ii + 1));
                    S2 = params.k2 * (1 - (q(ii + 1) / v(ii + 1)));
                    S3 = params.k3 * (1 - v(ii + 1));
                    S(ii + 1) = V0 * (S1 + S2 + S3);
                    OEF(ii + 1) = 1 - (1 - E0) .^ (1. / ii_flow);
                    CMRO2(ii + 1) = (OEF(ii + 1) / E0) * IN_FLOW(ii);
                    IN_FLOW(ii + 1) = ii_flow;
                end
                Sr = convolve_vecs(S(1:rd{rr, ss} / params.delta_t)', 1, fs, 1 / tr);
                run_preds{rr, ss}(:, pp) = Sr;
            end
        end
        session_run_preds = run_preds(:, ss); 
        save(fpath, 'session_run_preds', '-v7.3');
        fprintf('\n');
    end
end
run_preds(empty_cells) = {[]};
model.run_preds = run_preds;

end
