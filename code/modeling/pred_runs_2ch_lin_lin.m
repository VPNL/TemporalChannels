function model = pred_runs_2ch_lin_lin(model)
% Generates run predictors using the 2 temporal-channel with linear 
% sustained and linear transient channel. 

% get design parameters
fs = model.fs; tr = model.tr; stim = model.stim;
nruns_max = size(stim, 1); empty_cells = cellfun(@isempty, stim);
irfs_names = fieldnames(model.irfs); irfs = [];
for ff = 1:length(irfs_names)
    iname = model.irfs.(irfs_names{ff});
    irfs.(irfs_names{ff}) = repmat(iname, nruns_max, 1);
end

% generate run predictors for each session
predS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    stim, irfs.nrfS, 'uni', false); predS(empty_cells) = {[]};
predT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, fs), ...
    stim, irfs.nrfT, 'uni', false); predT(empty_cells) = {[]};
fmriS = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predS, irfs.hrf, 'uni', false); fmriS(empty_cells) = {[]};
fmriT = cellfun(@(X, Y) convolve_vecs(X, Y, fs, 1 / tr), ...
    predT, irfs.hrf, 'uni', false); fmriT(empty_cells) = {[]};
run_preds = cellfun(@(X, Y) [X Y * model.normT], fmriS, fmriT, 'uni', false);
model.run_preds = run_preds;

end
