function model = pred_runs_3ch_exp_cquad_rect(model)
% Generates run predictors using the 3 temporal-channel model with adapted
% sustained and quadratic + compressed transient and persistent channels.

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
predTq = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs) .^ 2, ...
    stim, irfs.nrfT, 'uni', false); predTq(empty_cells) = {[]};
predTc = cellfun(@(X) rectify(X, 'abs', .001) .^ .1, predTq, 'uni', false);
%predP = cellfun(@(X, Y) rectify(convolve_vecs(X, Y, fs, fs)) .^ 2, ...
%    stim, irfs.nrfP, 'uni', false); predP(empty_cells) = {[]};
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    model.adapt_act, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predTc, irfs.hrf, 'uni', false); fmriT(empty_cells) = {[]};
%fmriP = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
%    predP, irfs.hrf, 'uni', false); fmriP(empty_cells) = {[]};
fmriP = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    model.persist_act, irfs.hrf, 'uni', false); fmriP(empty_cells) = {[]};
model.run_preds = cellfun(@(X, Y, Z) [X Y * model.normT Z * model.normP], ...
    fmriS, fmriT, fmriP, 'uni', false);

end
