function obj_fun = tch_obj_fun_1ch_exp(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 1ch-exp model (single channel with exponential adaptation decay)
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

if ~strcmp(model.type, '1ch-exp'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;

adapt_fun = @(y) exp(-(1:60000) / y);
conv_sn = @(x, y) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@(XX) convolve_vecs(XX, irfs.nrfS{1}, 1, 1), x, 'uni', false), ...
    repmat({adapt_fun(y)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
conv_nb = @(x, y) cellfun(@(N) convolve_vecs(N, irfs.hrf{1}, fs, 1 / tr), ...
    conv_sn(x, y), 'uni', false);
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
comp_ws = @(x, y, m, b0) cell2mat(conv_nb(x, y)) \ cell2mat(comp_bs(m, b0));
pred_bs = @(x, y, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(x, y), repmat({comp_ws(x, y, m, b0)'}, nruns, 1), 'uni', false);
calc_br = @(x, y, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(x, y, m, b0), comp_bs(m, b0), 'uni', false);
calc_me = @(x, y, m, b0) sum(cell2mat(calc_br(x, y, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), run_avgs, baseline);

end
