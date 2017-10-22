function obj_fun = tch_obj_fun_1ch_exp(roi, model)

if ~strcmp(model.type, '1ch-exp'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); npreds = size(stim{1}, 2);
irfs = model.irfs; fs = model.fs; tr = model.tr;
run_avgs = roi.run_avgs; baseline = roi.baseline;
param_names = fieldnames(model.params); nparams = length(param_names);

adapt_fun = @(y) exp(-(1:60000) / y);
conv_sn = @(x, y) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@(XX) convolve_vecs(XX, irfs.nrfS{1}, 1, 1), x, 'uni', false), ...
    repmat({adapt_fun(y)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
conv_nb = @(x, y) cellfun(@(N) convolve_vecs(N, irfs.hrf{1}, fs, 1 / tr), ...
    conv_sn(x, y), 'uni', false);
pred_bs = @(x, y, b) cellfun(@(P, B) P .* repmat(B, size(P, 1), 1), ...
    conv_nb(x, y), repmat({b}, nruns, 1), 'uni', false);
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
calc_br = @(x, y, b, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, b), comp_bs(m, b0), 'uni', false);
calc_me = @(x, y, b, m, b0) sum(cell2mat(calc_br(x, y, b, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x([1:npreds] + nparams), run_avgs, baseline);

end
