function model = pred_runs_standard(model)
% Generates run predictors using the standard linear model. 

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
run_preds = cellfun(@(X) zeros(X/tr, ncats), rd, 'uni', false);
empty_cells = cellfun(@isempty, run_preds); run_preds(empty_cells) = {[]};
run_preds = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), stim, irfs.hrf, 'uni', false);
model.run_preds = run_preds;

end
