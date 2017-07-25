function model = run_preds_all(model)

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
if strcmp(model.type, 'standard')
    % if using standard model
    pred_mat = cellfun(@(X, Y) convolve_vecs(X, fs, tr, Y), stim, irfs.hrf, 'uni', false);
elseif sum(strcmp(model.type, {'2ch' '2ch-cts' '2ch-dcts' '3ch'}))
    % if using a multi-channel model
    predS = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y), stim, irfs.nrfS, 'uni', false);
    predT = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y).^2, stim, irfs.nrfT, 'uni', false);
    switch model.type
        case '2ch-cts'
            predS = cellfun(@(X, Y) X.^Y, predS, params.e, 'uni', false);
        case '2ch-dcts'
            predSn = cellfun(@(X) X.^2, predS, 'uni', false);
            predSd = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y), predS, irfs.lpf, 'uni', false);
            predSd = cellfun(@(X, Y) X.^2 + Y.^2, predSd, params.sigma, 'uni', false);
            predS = cellfun(@(X, Y) X./Y, predSn, predSd, 'uni', false);
        case '3ch'
            csD = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y), stimD, irfs.nrfD, 'uni', false);
            predDn = cellfun(@(X) X.^2, csD, 'uni', false);
            predDd = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y), csD, irfs.lpf, 'uni', false);
            predDd = cellfun(@(X, Y) X.^2 + Y.^2, predDd, params.sigma, 'uni', false);
            predD = cellfun(@(X, Y) X./Y, predDn, predDd, 'uni', false);
            predD(empty_cells) = {1};
            fmriD = cellfun(@(X, Y) convolve_vecs(X, fs, tr, Y), predD, irfs.hrf, 'uni', false);
            fmriD(empty_cells) = {[]};
            predS = cellfun(@(X, Y) X.^Y, predS, params.e, 'uni', false);
    end
    predS(empty_cells) = {1}; predT(empty_cells) = {1};
    % convolve neural predictors with HRF
    fmriS = cellfun(@(X, Y) convolve_vecs(X, fs, tr, Y), predS, irfs.hrf, 'uni', false);
    fmriT = cellfun(@(X, Y) convolve_vecs(X, fs, tr, Y), predT, irfs.hrf, 'uni', false);
    fmriS(empty_cells) = {[]}; fmriT(empty_cells) = {[]};
    % order fMRI predictors in the same design matrix
    if strcmp(model.type, '3ch')
        pred_mat = cellfun(@(X, Y, Z) [X Y*model.normT Z*model.normD], fmriS, fmriT, fmriD, 'uni', false);
    else
        pred_mat = cellfun(@(X, Y) [X Y*model.normT], fmriS, fmriT, 'uni', false);
    end
elseif strcmp(model.type, 'cts')
    % if using standard CTS model
    predS = cellfun(@(X, Y, Z) convolve_vecs(X, fs, 1/fs, Y).^Z, stim, irfs.nrfS, params.e, 'uni', false);
    predS(empty_cells) = {[]};
    pred_mat = cellfun(@(X, Y) convolve_vecs(X, model.fs, tr, Y), predS, irfs.hrf, 'uni', false);
elseif strcmp(model.type, 'dcts')
    % if using standard dCTS model
    predS = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y).^2, stim, irfs.nrfS, 'uni', false);
    predSn = cellfun(@(X) X.^2, predS, 'uni', false);
    predSd = cellfun(@(X, Y) convolve_vecs(X, fs, 1/fs, Y), predS, irfs.lpf, 'uni', false);
    predSd = cellfun(@(X, Y) X.^2 + Y.^2, predSd, params.sigma, 'uni', false);
    predS = cellfun(@(X, Y) X./Y, predSn, predSd, 'uni', false);
    predS(empty_cells) = {[]};
    pred_mat = cellfun(@(X, Y) convolve_vecs(X, fs, tr, Y), predS, irfs.hrf, 'uni', false);
end

model.pred_mat = pred_mat;

end
