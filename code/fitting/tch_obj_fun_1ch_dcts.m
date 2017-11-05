function obj_fun = tch_obj_fun_1ch_dcts(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 1ch-dcts model (single channel with dynamtic CTS)
% 
% INPUTS:
%   1) roi: tchROI object containing single session
%   2) model: tchModel object for the same session
% 
% OUTPUTS:
%   obj_fun: anonymous objective function in the form of y = f(x0), where
%   x0 is a vector of parameters to evaluate and y is the sum of squared
%   residual error between model predictions and run response time series
% 
% AS 10/2017

if ~strcmp(model.type, '1ch-dcts'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;

irf_fun = @(x) (0:999) .* exp(-(0:999) / (x * 1000));
div_lpf = @(f) exp(-(0:999) / f) / sum(exp(-(0:999) / f));
conv_sn = @(x, y) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    x, repmat({irf_fun(y)}, nruns, 1), 'uni', false);
conv_nf = @(x, y, f) cellfun(@(X, F) convolve_vecs(X, F, fs, fs) .^ 2, ...
    conv_sn(x, y), repmat({div_lpf(f)}, nruns, 1), 'uni', false);
comp_dn = @(x, y, z, f) cellfun(@(N, F, Z) (N .^ 2) ./ (F + Z .^ 2), ...
    conv_sn(x, y), conv_nf(x, y, f), repmat({z}, nruns, 1), 'uni', false);
conv_nb = @(x, y, z, f) cellfun(@(N) convolve_vecs(N, irfs.hrf{1}, fs, 1 / tr), ...
    comp_dn(x, y, z, f), 'uni', false);
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
comp_ws = @(x, y, z, f, m, b0) cell2mat(conv_nb(x, y, z, f)) \ cell2mat(comp_bs(m, b0));
pred_bs = @(x, y, z, f, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(x, y, z, f), repmat({comp_ws(x, y, z, f, m, b0)'}, nruns, 1), 'uni', false);
calc_br = @(x, y, z, f, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, z, f, m, b0), comp_bs(m, b0), 'uni', false);
calc_me = @(x, y, z, f, m, b0) sum(cell2mat(calc_br(x, y, z, f, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), run_avgs, baseline);

end
