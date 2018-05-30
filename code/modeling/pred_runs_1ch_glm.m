function model = pred_runs_1ch_glm(model)
% Generates run predictors using a simple general linear model. 

% get design parameters
fs = model.fs; tr = model.tr; stim = model.stim;
nruns_max = size(stim, 1); empty_cells = cellfun(@isempty, stim);
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end

% generate run predictors for each session
stim_trs = cellfun(@(X) reshape(X', size(X, 2), tr * fs, []), ...
    stim, 'uni', false); stim_trs(empty_cells) = {[]};
stim_trs = cellfun(@(X) ceil(squeeze(mean(X, 2))'), ...
    stim_trs, 'uni', false); stim_trs(empty_cells) = {[]};
run_preds = cellfun(@(X) convolve_vecs(X, canonical_hrf(tr), 1, 1), ...
    stim_trs, 'uni', false); run_preds(empty_cells) = {[]};
model.run_preds = run_preds;

end
