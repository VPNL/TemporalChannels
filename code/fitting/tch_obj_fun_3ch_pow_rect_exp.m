function obj_fun = tch_obj_fun_3ch_pow_rect_exp(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 3ch-pow-rect-exp model (3-channel model with CTS-pow on sustained,
% rectified transient, and exponential persistent channels).
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

if ~strcmp(model.type, '3ch-pow-rect-exp'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
persist_fun = @(z) exp(-(1:12000) / (z * 1000));
% sustained response: (stimulus * sustained IRF)^epsilon
conv_snS = @(x, y) cellfun(@(X, Y) convolve_vecs(X, irfs.nrfS{1}, 1, 1) .^ Y, ...
    x, repmat({y}, nruns, 1), 'uni', false);
% transient response: max(0, stimulus * transient IRF)
conv_snT = @(x) cellfun(@(X) rectify(convolve_vecs(X, irfs.nrfT{1}, 1, 1), 'positive'), ...
    x, 'uni', false);
% persistent response: persistent function x exponential[tau_pe]
poffsets = cellfun(@(X, Y) [X(2:end) Y], model.onsets, model.run_durs, 'uni', false);
conv_snP = @(x, z) cellfun(@(X, Z, ON, OFF) code_exp_decay(X, ON, OFF, Z, fs), ...
    cellfun(@code_persist_act, x, 'uni', false), repmat({persist_fun(z)}, nruns, 1), ...
    model.offsets, poffsets, 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(x, y) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(x, y), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(x) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(x), 'uni', false);
% persistent BOLD: persistent response * HRF
conv_nbP = @(x, z) cellfun(@(NP) convolve_vecs(NP, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snP(x, z), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD, persistent BOLD]
conv_nb = @(x, y, z) cellfun(@(S, T, P) [S T P], ...
    conv_nbS(x, y), conv_nbT(x), conv_nbP(x, z), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(x, y, z, m, b0) cell2mat(conv_nb(x, y, z)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(x, y, z, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(x, y, z), repmat({comp_ws(x, y, z, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(x, y, z, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, z, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(x, y, z, m, b0) sum(cell2mat(calc_br(x, y, z, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), run_avgs, baseline);

end
