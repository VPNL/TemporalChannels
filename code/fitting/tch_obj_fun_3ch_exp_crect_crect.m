function obj_fun = tch_obj_fun_3ch_exp_crect_crect(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% 2-channel model with optimized adapated sustained and separate rectifed + 
% compressed transient channels for on and off responses).
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

if ~strcmp(model.type, '3ch-exp-crect-crect'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
nrfS_fun = @(tau_s) tch_irfs('S', tau_s);
nrfT_fun = @(tau_t) tch_irfs('T', tau_t);
adapt_fun = @(tau_ae) exp(-(1:60000) / (tau_ae * 10000));
% sustained response: (stimulus * sustained IRF) x exponential[tau_ae]
conv_snS = @(s, tau_s, tau_ae) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@(XX, YY) convolve_vecs(XX, YY, 1, 1), s, repmat({nrfS_fun(tau_s)}, nruns, 1), 'uni', false), ...
    repmat({adapt_fun(tau_ae)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
% transient on response: max(0, (stimulus * transient IRF))^2)^epsilon
conv_snTr1 = @(s, tau_t) cellfun(@(X, Y) rectify(convolve_vecs(X, Y, 1, 1), 'positive'), ...
    s, repmat({nrfT_fun(tau_t)}, nruns, 1), 'uni', false);
conv_snTc1 = @(s, tau_t) cellfun(@(X, Y) rectify(X .^ 2, 'abs', .001) .^ .1, ...
    conv_snTr1(s, tau_t), 'uni', false);
% transient off response: min(0, (stimulus * transient IRF))^2)^epsilon
conv_snTr2 = @(s, tau_t) cellfun(@(X, Y) rectify(convolve_vecs(X, Y, 1, 1), 'negative'), ...
    s, repmat({nrfT_fun(tau_t)}, nruns, 1), 'uni', false);
conv_snTc2 = @(s, tau_t) cellfun(@(X, Y) rectify(X .^ 2, 'abs', .001) .^ .1, ...
    conv_snTr2(s, tau_t), 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(s, tau_s, tau_ae) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(s, tau_s, tau_ae), 'uni', false);
% transient on BOLD: transient response * HRF
conv_nbT1 = @(s, tau_t) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snTc1(s, tau_t), 'uni', false);
% transient off BOLD: transient response * HRF
conv_nbT2 = @(s, tau_t) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snTc2(s, tau_t), 'uni', false);
% channel predictors: [sustained BOLD, transient on BOLD, transient off BOLD]
conv_nb = @(s, tau_s, tau_t, tau_ae) cellfun(@(S, T1, T2) [S T1 T2], ...
    conv_nbS(s, tau_s, tau_ae), conv_nbT1(s, tau_t), conv_nbT2(s, tau_t), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(s, tau_s, tau_t, tau_ae, m, b0) ...
    cell2mat(conv_nb(s, tau_s, tau_t, tau_ae)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(s, tau_s, tau_t, tau_ae, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(s, tau_s, tau_t, tau_ae), repmat({comp_ws(s, tau_s, tau_t, tau_ae, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(s, tau_s, tau_t, tau_ae, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, tau_s, tau_t, tau_ae, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(s, tau_s, tau_t, tau_ae, m, b0) ...
    sum(cell2mat(calc_br(s, tau_s, tau_t, tau_ae, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), run_avgs, baseline);

end
