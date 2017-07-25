function model = run_preds_3ch(model)
% Generates run predictors using the 3 temporal-channel model. 

% get design parameters
params_init = model.params; irfs_init = model.irfs;
fs = model.fs; tr = model.tr; rd = model.run_durs;
stim = model.stim;  stimD = model.stimD;
cat_list = unique([model.cats{:}]); ncats = length(cat_list);
[nruns_max, ~] = size(model.stimfiles);
params_names = fieldnames(params_init); params = [];
for pp = 1:length(params_names)
    params.(params_names{pp}) = repmat(params_init.(params_names{pp}), nruns_max, 1);
end
irfs_names = fieldnames(irfs_init); irfs = [];
for ff = 1:length(irfs_names)
    irfs.(irfs_names{ff}) = repmat(irfs_init.(irfs_names{ff}), nruns_max, 1);
end

% generate run predictors for each session
pred_mat = cellfun(@(X) zeros(X/model.tr, ncats), rd, 'uni', false);
empty_cells = cellfun(@isempty, pred_mat); pred_mat(empty_cells) = {[]};
predS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), stim, irfs.nrfS, 'uni', false);
predT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs) .^ 2, stim, irfs.nrfT, 'uni', false);
predS(empty_cells) = {1}; predT(empty_cells) = {1};
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predS, irfs.hrf, 'uni', false);
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predT, irfs.hrf, 'uni', false);
fmriS(empty_cells) = {[]}; fmriT(empty_cells) = {[]};
% code delay channel
csD = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), stimD, irfs.nrfD, 'uni', false);
predDn = cellfun(@(X) X .^ 2, csD, 'uni', false);
predDd = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), csD, irfs.lpf, 'uni', false);
predDd = cellfun(@(X, Y) X .^ 2 + Y .^ 2, predDd, params.sigma, 'uni', false);
predD = cellfun(@(X, Y) X ./ Y, predDn, predDd, 'uni', false);
predD(empty_cells) = {1};
fmriD = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predD, irfs.hrf, 'uni', false);
fmriD(empty_cells) = {[]};
pred_mat = cellfun(@(X, Y, Z) [X Y * model.normT Z * model.normD], fmriS, fmriT, fmriD, 'uni', false);
model.pred_mat = pred_mat;

end
