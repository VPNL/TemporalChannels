function obj_fun = tch_obj_fun_2ch_lin_quad(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 2ch-lin-quad model (optimized 2-channel model with sustained and  
% quadratic transient channels).
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

if ~strcmp(model.type, '2ch-lin-quad'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
nrfS_fun = @(tau_s) tch_irfs('S', tau_s);
nrfT_fun = @(tau_s) tch_irfs('T', tau_s);
% sustained response: stimulus * sustained IRF
conv_snS = @(s, tau_s) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    s, repmat({nrfS_fun(tau_s)}, nruns, 1), 'uni', false);
% transient response: (stimulus * transient IRF)^2
conv_snT = @(s, tau_s) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1) .^ 2, ...
    s, repmat({nrfT_fun(tau_s)}, nruns, 1), 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(s, tau_s) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(s, tau_s), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(s, tau_s) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(s, tau_s), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD, persistent BOLD]
conv_nb = @(s, tau_s) cellfun(@(S, T) [S T], ...
    conv_nbS(s, tau_s), conv_nbT(s, tau_s), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(s, tau_s, m, b0) cell2mat(conv_nb(s, tau_s)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(s, tau_s, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(s, tau_s), repmat({comp_ws(s, tau_s, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(s, tau_s, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, tau_s, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(s, tau_s, m, b0) sum(cell2mat(calc_br(s, tau_s, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), run_avgs, baseline);

end