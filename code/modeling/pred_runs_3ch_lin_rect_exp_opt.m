function model = pred_runs_3ch_lin_rect_exp_opt(model)
% Generates run predictors using the 3 temporal-channel model with 
% optimized linear sustained, rectified transient, and delay channel.

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
    stim, irfs.nrfS, 'uni', false); predS(empty_cells) = {1};
predTr = cellfun(@(X, Y) rectify(convolve_vecs(X, Y, fs, fs)), ...
    stim, irfs.nrfT, 'uni', false); predTr(empty_cells) = {1};
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predS, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predTr, irfs.hrf, 'uni', false); fmriT(empty_cells) = {[]};
fmriD = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    model.delay_act, irfs.hrf, 'uni', false); fmriD(empty_cells) = {[]};
run_preds = cellfun(@(X, Y, Z) [X Y * model.normT Z * model.normD], ...
    fmriS, fmriT, fmriD, 'uni', false); run_preds(empty_cells) = {[]};
model.run_preds = run_preds;

end
