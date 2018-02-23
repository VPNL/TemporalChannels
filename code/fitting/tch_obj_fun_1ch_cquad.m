function obj_fun = tch_obj_fun_1ch_cquad(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 1ch-cquad model (optimized 1-channel model with quadratic +
% compressed transient channels).
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
% AS 1/2018

if ~strcmp(model.type, '1ch-cquad'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
nrfT_fun = @(tau_t) tch_irfs('T', tau_t);
% transient response: ((stimulus * transient IRF)^2)^epsilon
conv_snTq = @(s, tau_t) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1) .^ 2, ...
    s, repmat({nrfT_fun(tau_t)}, nruns, 1), 'uni', false);
conv_snTc = @(s, tau_t) cellfun(@(X, Y) rectify(X, 'abs', .001) .^ .1, ...
    conv_snTq(s, tau_t), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(s, tau_t) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snTc(s, tau_t), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD, persistent BOLD]
conv_nb = @(s, tau_t) cellfun(@(T) [T], ...
    conv_nbT(s, tau_t), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(s, tau_t, m, b0) cell2mat(conv_nb(s, tau_t)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(s, tau_t, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(s, tau_t), repmat({comp_ws(s, tau_t, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(s, tau_t, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, tau_t, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(s, tau_t, m, b0) sum(cell2mat(calc_br(s, tau_t, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), run_avgs, baseline);

end
