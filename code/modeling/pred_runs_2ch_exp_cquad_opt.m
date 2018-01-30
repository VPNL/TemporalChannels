function model = pred_runs_2ch_exp_cquad_opt(model)
% Generates run predictors using the 2ch-exp-cquad-opt model (2-channel
% model with optimized adapted sustained and compressed quadratic transient
% channels).

% get design parameters
fs = model.fs; tr = model.tr; stim = model.stim;
nruns_max = size(stim, 1); empty_cells = cellfun(@isempty, stim);
params_names = fieldnames(model.params); params = [];
for pp = 1:length(params_names)
    pname = model.params.(params_names{pp});
    params.(params_names{pp}) = repmat(pname, nruns_max, 1);
end
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end

% generate run predictors for each session
predS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    model.adapt_act, irfs.nrfS, 'uni', false); predS(empty_cells) = {1};
predTl = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    stim, irfs.nrfT, 'uni', false);
predTn = cellfun(@(X) X .^ 2, predTl, 'uni', false);
predTf = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    predTl, irfs.lpf, 'uni', false);
predTd = cellfun(@(X, Y) X .^ 2 + Y .^ 2, ...
    predTf, params.sigma, 'uni', false);
predTcq = cellfun(@(X, Y) X ./ Y, predTn, predTd, 'uni', false);
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predS, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predTcq, irfs.hrf, 'uni', false); fmriT(empty_cells) = {[]};
run_preds = cellfun(@(X, Y) [X Y * model.normT], ...
    fmriS, fmriT, 'uni', false); run_preds(empty_cells) = {[]};
model.run_preds = run_preds;

end
