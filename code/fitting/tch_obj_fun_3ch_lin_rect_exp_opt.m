function obj_fun = tch_obj_fun_3ch_lin_rect_exp_opt(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 3ch-lin-quad-exp model (3-channel model with optimized linear 
% sustained, rectified transient, and exponential persistent channels).
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

if ~strcmp(model.type, '3ch-lin-rect-exp-opt'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
nrfS_fun = @(tau, n1, n2, k) tch_irfs('S', tau, n1, n2, k, fs);
nrfT_fun = @(tau, n1, n2, k) tch_irfs('T', tau, n1, n2, k, fs);
persist_fun = @(tau_pe) exp(-(1:12000) / (tau_pe * 1000));
% sustained response: stimulus * sustained IRF
conv_snS = @(s, tau, n1, n2, k) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    s, repmat({nrfS_fun(tau, n1, n2, k)}, nruns, 1), 'uni', false);
% transient response: max(0, stimulus * transient IRF)
conv_snT = @(s, tau, n1, n2, k) cellfun(@(X, Y) rectify(convolve_vecs(X, Y, 1, 1)), ...
    s, repmat({nrfT_fun(tau, n1, n2, k)}, nruns, 1), 'uni', false);
% persistent response: persistent function x exponential[tau_pe]
poffsets = cellfun(@(X, Y) [X(2:end) Y], model.onsets, model.run_durs, 'uni', false);
conv_snP = @(s, tau_pe) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@code_persist_act, s, 'uni', false), repmat({persist_fun(tau_pe)}, nruns, 1), ...
    model.offsets, poffsets, 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(s, tau, n1, n2, k) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(s, tau, n1, n2, k), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(s, tau, n1, n2, k) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(s, tau, n1, n2, k), 'uni', false);
% persistent BOLD: persistent response * HRF
conv_nbP = @(s, tau_pe) cellfun(@(NP) convolve_vecs(NP, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snP(s, tau_pe), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD, persistent BOLD]
conv_nb = @(s, tau, n1, n2, k, tau_pe) cellfun(@(S, T, P) [S T P], ...
    conv_nbS(s, tau, n1, n2, k), conv_nbT(s, tau, n1, n2, k), conv_nbP(s, tau_pe), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(s, tau, n1, n2, k, tau_pe, m, b0) ...
    cell2mat(conv_nb(s, tau, n1, n2, k, tau_pe)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(s, tau, n1, n2, k, tau_pe, m, b0) ...
    cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), conv_nb(s, tau, n1, n2, k, tau_pe), ...
    repmat({comp_ws(s, tau, n1, n2, k, tau_pe, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(s, tau, n1, n2, k, tau_pe, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, tau, n1, n2, k, tau_pe, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(s, tau, n1, n2, k, tau_pe, m, b0) sum(cell2mat(calc_br(s, tau, n1, n2, k, tau_pe, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), x(4), x(5), run_avgs, baseline);

end
