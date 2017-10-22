function obj_fun = tch_obj_fun_1ch_pow(roi, model)

if ~strcmp(model.type, '1ch-pow'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); npreds = size(stim{1}, 2);
irfs = model.irfs; fs = model.fs; tr = model.tr;
run_avgs = roi.run_avgs; baseline = roi.baseline;
param_names = fieldnames(model.params); nparams = length(param_names);

irf_fun = @(x) (0:999) .* exp(-(0:999) / x);
conv_sn = @(x, y, z) cellfun(@(X, Y, Z) convolve_vecs(X, Y, 1, 1) .^ Z, ...
    x, repmat({irf_fun(y)}, nruns, 1), repmat({z}, nruns, 1), 'uni', false);
conv_nb = @(x, y, z) cellfun(@(N) convolve_vecs(N, irfs.hrf{1}, fs, 1 / tr), ...
    conv_sn(x, y, z), 'uni', false);
pred_bs = @(x, y, z, b) cellfun(@(P, B) P .* repmat(B, size(P, 1), 1), ...
    conv_nb(x, y, z), repmat({b}, nruns, 1), 'uni', false);
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
calc_br = @(x, y, z, b, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, z, b), comp_bs(m, b0), 'uni', false);
calc_me = @(x, y, z, b, m, b0) sum(cell2mat(calc_br(x, y, z, b, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x([1:npreds] + nparams), run_avgs, baseline);

end
