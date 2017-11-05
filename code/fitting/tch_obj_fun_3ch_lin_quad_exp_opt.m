function obj_fun = tch_obj_fun_3ch_lin_quad_exp_opt(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 3ch-lin-quad-exp model (3-channel model with optimized linear 
% sustained, quadratic transient, and exponential delay channels).
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

if ~strcmp(model.type, '3ch-lin-quad-exp-opt'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
nrfS_fun = @(x, y) tch_irfs('S', x, y, fs);
nrfT_fun = @(x, y) tch_irfs('T', x, y, fs);
delay_fun = @(z) exp(-(1:12000) / (z * 1000));
% sustained response: stimulus * sustained IRF
conv_snS = @(s, x, y) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    s, repmat({nrfS_fun(x, y)}, nruns, 1), 'uni', false);
% transient response: (stimulus * transient IRF)^2
conv_snT = @(s, x, y) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1) .^ 2, ...
    s, repmat({nrfT_fun(x, y)}, nruns, 1), 'uni', false);
% delay response: delay function x exponential[tau_de]
doffsets = cellfun(@(X, Y) [X(2:end) Y], model.onsets, model.run_durs, 'uni', false);
conv_snD = @(s, z) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@code_delay_act, s, 'uni', false), repmat({delay_fun(z)}, nruns, 1), ...
    model.offsets, doffsets, 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(s, x, y) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(s, x, y), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(s, x, y) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(s, x, y), 'uni', false);
% delay BOLD: delay response * HRF
conv_nbD = @(s, z) cellfun(@(ND) convolve_vecs(ND, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snD(s, z), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD, delay BOLD]
conv_nb = @(s, x, y, z) cellfun(@(S, T, D) [S T D], ...
    conv_nbS(s, x, y), conv_nbT(s, x, y), conv_nbD(s, z), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(s, x, y, z, m, b0) cell2mat(conv_nb(s, x, y, z)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(s, x, y, z, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(s, x, y, z), repmat({comp_ws(s, x, y, z, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(s, x, y, z, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, x, y, z, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(s, x, y, z, m, b0) sum(cell2mat(calc_br(s, x, y, z, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), run_avgs, baseline);

end
