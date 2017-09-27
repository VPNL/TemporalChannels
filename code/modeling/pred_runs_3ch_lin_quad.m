function model = pred_runs_3ch_lin_quad(model)
% Generates run predictors using the 3 temporal-channel model with linear
% sustained, quradratic transient, and optimized delay channel.

% get design parameters
params_init = model.params; irfs_init = model.irfs;
fs = model.fs; tr = model.tr; rd = model.run_durs; stim = model.stim;
cat_list = unique([model.cats{:}]); ncats = length(cat_list);
[nruns_max, ~] = size(model.run_durs);
params_names = fieldnames(params_init); params = [];
for pp = 1:length(params_names)
    params.(params_names{pp}) = repmat(params_init.(params_names{pp}), nruns_max, 1);
end
irfs_names = fieldnames(irfs_init); irfs = [];
for ff = 1:length(irfs_names)
    irfs.(irfs_names{ff}) = repmat(irfs_init.(irfs_names{ff}), nruns_max, 1);
end

% generate run predictors for each session
run_preds = cellfun(@(X) zeros(X / tr, ncats), rd, 'uni', false);
empty_cells = cellfun(@isempty, run_preds); run_preds(empty_cells) = {[]};
predS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), stim, irfs.nrfS, 'uni', false);
predT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), stim, irfs.nrfT, 'uni', false);
predTq = cellfun(@(X) X .^ 2, predT, 'uni', false);
predD = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), stim, irfs.nrfD, 'uni', false);
predDr = cellfun(@(X) rectify(X, 'negative') .^ 2, predD, 'uni', false);
predS(empty_cells) = {1}; predTq(empty_cells) = {1}; predDr(empty_cells) = {1};
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predS, irfs.hrf, 'uni', false);
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predTq, irfs.hrf, 'uni', false);
fmriD = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predDr, irfs.hrf, 'uni', false);
fmriS(empty_cells) = {[]}; fmriT(empty_cells) = {[]}; fmriD(empty_cells) = {[]};
run_preds = cellfun(@(X, Y, Z) [X Y * model.normT Z * model.normD], fmriS, fmriT, fmriD, 'uni', false);
model.run_preds = run_preds;

end
