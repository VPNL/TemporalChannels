function model = pred_runs_3ch_exp_crect_crect(model)
% Generates run predictors using the 3 temporal-channel model with adapted
% sustained and separate rectified + compressed transient channels for on
% and off responses.

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
predTr1 = cellfun(@(X, Y) rectify(convolve_vecs(X, Y, fs, fs), 'positive') .^ 2, ...
    stim, irfs.nrfT, 'uni', false); predTr1(empty_cells) = {[]};
predTc1 = cellfun(@(X) rectify(X, 'abs', .001) .^ .1, predTr1, 'uni', false);
predTr2 = cellfun(@(X, Y) rectify(convolve_vecs(X, Y, fs, fs), 'negative') .^ 2, ...
    stim, irfs.nrfT, 'uni', false); predTr2(empty_cells) = {[]};
predTc2 = cellfun(@(X) rectify(X, 'abs', .001) .^ .1, predTr2, 'uni', false);
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    model.adapt_act, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
fmriT1 = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predTc1, irfs.hrf, 'uni', false); fmriT1(empty_cells) = {[]};
fmriT2 = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predTc2, irfs.hrf, 'uni', false); fmriT2(empty_cells) = {[]};
model.run_preds = cellfun(@(X, Y, Z) [X Y * model.normT Z * model.normT], ...
    fmriS, fmriT1, fmriT2, 'uni', false);

end
