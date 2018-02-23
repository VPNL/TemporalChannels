function obj_fun = tch_obj_fun_2ch_exp_sig(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% 2-channel model with optimized adapated sustained and sigmoid transient 
% channels).
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

if ~strcmp(model.type, '2ch-exp-sig'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
nrfS_fun = @(tau_s) tch_irfs('S', tau_s);
nrfT_fun = @(tau_s) tch_irfs('T', tau_s);
adapt_fun = @(tau_ae) exp(-(1:60000) / (tau_ae * 10000));
% sustained response: (stimulus * sustained IRF) x exponential[tau_ae]
conv_snS = @(s, tau_s, tau_ae) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@(XX, YY) convolve_vecs(XX, YY, 1, 1), s, repmat({nrfS_fun(tau_s)}, nruns, 1), 'uni', false), ...
    repmat({adapt_fun(tau_ae)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
% transient response: tch_sigmoid(stimulus * transient IRF)
conv_snT = @(s, tau_s) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    s, repmat({nrfT_fun(tau_s)}, nruns, 1), 'uni', false);
conv_snTs = @(s, tau_s, Lp, Kp, Ln, Kn) cellfun(@(X, lp, kp, ln, kn) tch_sigmoid(X, lp, kp, ln, kn), ...
    conv_snT(s, tau_s), repmat({Lp}, nruns, 1), repmat({Kp}, nruns, 1), repmat({Ln}, nruns, 1), repmat({Kn}, nruns, 1), 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(s, tau_s, tau_ae) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(s, tau_s, tau_ae), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(s, tau_s, Lp, Kp, Ln, Kn) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snTs(s, tau_s, Lp, Kp, Ln, Kn), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD]
conv_nb = @(s, tau_s, tau_ae, Lp, Kp, Ln, Kn) cellfun(@(S, T) [S T], ...
    conv_nbS(s, tau_s, tau_ae), conv_nbT(s, tau_s, Lp, Kp, Ln, Kn), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: channel predictors \ measured signal
comp_ws = @(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0) ...
    cell2mat(conv_nb(s, tau_s, tau_ae, Lp, Kp, Ln, Kn)) \ cell2mat(comp_bs(m, b0));
% predicted signal: channel predictors x channel weights
pred_bs = @(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(s, tau_s, tau_ae, Lp, Kp, Ln, Kn), repmat({comp_ws(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0) ...
    sum(cell2mat(calc_br(s, tau_s, tau_ae, Lp, Kp, Ln, Kn, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), x(4), x(5), x(6), run_avgs, baseline);

end
