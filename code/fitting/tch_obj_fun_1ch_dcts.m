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
% generate IRFs/filters for optimization
irf_fun = @(x) (0:999) .* exp(-(0:999) / (x * 1000));
div_lpf = @(f) exp(-(0:1999) / f) / sum(exp(-(0:999) / f));
% linear response: stim * IRF[tau1]
conv_sn = @(x, y) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    x, repmat({irf_fun(y)}, nruns, 1), 'uni', false);
% filtered response: (linear response * low-pass filter[tau2])^2
conv_nf = @(x, y, f) cellfun(@(X, F) convolve_vecs(X, F, fs, fs) .^ 2, ...
    conv_sn(x, y), repmat({div_lpf(f)}, nruns, 1), 'uni', false);
% neural response: (linear response)^2 / (sigma^2 + filtered response)
comp_dn = @(x, y, z, f) cellfun(@(N, F, Z) (N .^ 2) ./ (F + Z .^ 2), ...
    conv_sn(x, y), conv_nf(x, y, f), repmat({z}, nruns, 1), 'uni', false);
% bold response: neural response * HRF
conv_nb = @(x, y, z, f) cellfun(@(N) convolve_vecs(N, irfs.hrf{1}, fs, 1 / tr), ...
    comp_dn(x, y, z, f), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: bold response \ measured signal
comp_ws = @(x, y, z, f, m, b0) cell2mat(conv_nb(x, y, z, f)) \ cell2mat(comp_bs(m, b0));
% predicted signal: bold response x channel weights
pred_bs = @(x, y, z, f, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(x, y, z, f), repmat({comp_ws(x, y, z, f, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(x, y, z, f, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, z, f, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(x, y, z, f, m, b0) sum(cell2mat(calc_br(x, y, z, f, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), run_avgs, baseline);

end
