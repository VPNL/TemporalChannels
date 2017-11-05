function obj_fun = tch_obj_fun_2ch_exp_rect(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 2ch-exp-rect model (2-channel model with adapted sustained and
% rectified transient channels).
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

if ~strcmp(model.type, '2ch-exp-rect'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;
% generate IRFs/filters for optimization
adapt_fun = @(y) exp(-(1:60000) / (y * 1000));
% sustained response: (stimulus * sustained IRF) x exponential[tau_ae]
conv_snS = @(x, y) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@(XX) convolve_vecs(XX, irfs.nrfS{1}, 1, 1), x, 'uni', false), ...
    repmat({adapt_fun(y)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
% transient response: max(0, stimulus * transient IRF)
conv_snT = @(x) cellfun(@(X) rectify(convolve_vecs(X, irfs.nrfT{1}, 1, 1), 'positive'), ...
    x, 'uni', false);
% sustained BOLD: sustained response * HRF
conv_nbS = @(x, y) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(x, y), 'uni', false);
% transient BOLD: transient response * HRF
conv_nbT = @(x) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(x), 'uni', false);
% channel predictors: [sustained BOLD, transient BOLD]
conv_nb = @(x, y) cellfun(@(S, T) [S T], conv_nbS(x, y), conv_nbT(x), 'uni', false);
% measured signal: time series - baseline estimates
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
% channel weights: bold response \ measured signal
comp_ws = @(x, y, m, b0) cell2mat(conv_nb(x, y)) \ cell2mat(comp_bs(m, b0));
% predicted signal: bold response x channel weights
pred_bs = @(x, y, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(x, y), repmat({comp_ws(x, y, m, b0)'}, nruns, 1), 'uni', false);
% model residuals: (predicted signal - measured signal)^2
calc_br = @(x, y, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, m, b0), comp_bs(m, b0), 'uni', false);
% model error: summed squared residuals for all run time series
calc_me = @(x, y, m, b0) sum(cell2mat(calc_br(x, y, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), run_avgs, baseline);

end
