function model = run_preds_cts(model)
% Generates run predictors using the CTS model. 

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
pred_mat = cellfun(@(X) zeros(X / model.tr, ncats), rd, 'uni', false);
empty_cells = cellfun(@isempty, pred_mat); pred_mat(empty_cells) = {[]};
predS = cellfun(@(X, Y, Z) convolve_vecs(X, Y, fs, fs) .^ Z, stim, irfs.nrfS, params.e, 'uni', false);
predS(empty_cells) = {[]};
pred_mat = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), predS, irfs.hrf, 'uni', false);
model.pred_mat = pred_mat;

end
