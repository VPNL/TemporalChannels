function obj_fun = tch_obj_fun_3ch_lin_rect_exp_opt(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 3ch-lin-quad-exp model (3-channel model with optimized linear 
% sustained, rectified transient, and exponential delay channels).
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

stim = model.stim; nruns = size(stim, 1); irfs = model.irfs; fs = model.fs;
run_avgs = roi.run_avgs; baseline = roi.baseline; tr = roi.tr;

nrfS_fun = @(x, y) tch_irfs('S', x, y, fs);
nrfT_fun = @(x, y) tch_irfs('T', x, y, fs);
delay_fun = @(z) exp(-(1:12000) / (z * 1000));
conv_snS = @(s, x, y) cellfun(@(X, Y) convolve_vecs(X, Y, 1, 1), ...
    s, repmat({nrfS_fun(x, y)}, nruns, 1), 'uni', false);
conv_snT = @(s, x, y) cellfun(@(X, Y) rectify(convolve_vecs(X, Y, 1, 1)), ...
    s, repmat({nrfT_fun(x, y)}, nruns, 1), 'uni', false);
doffsets = cellfun(@(X, Y) [X(2:end) Y], model.onsets, model.run_durs, 'uni', false);
conv_snD = @(s, z) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@code_delay_act, s, 'uni', false), repmat({delay_fun(z)}, nruns, 1), ...
    model.offsets, doffsets, 'uni', false);
conv_nbS = @(s, x, y) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(s, x, y), 'uni', false);
conv_nbT = @(s, x, y) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(s, x, y), 'uni', false);
conv_nbD = @(s, z) cellfun(@(ND) convolve_vecs(ND, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snD(s, z), 'uni', false);
conv_nb = @(s, x, y, z) cellfun(@(S, T, D) [S T D], ...
    conv_nbS(s, x, y), conv_nbT(s, x, y), conv_nbD(s, z), 'uni', false);
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
comp_ws = @(s, x, y, z, m, b0) cell2mat(conv_nb(s, x, y, z)) \ cell2mat(comp_bs(m, b0));
pred_bs = @(s, x, y, z, m, b0) cellfun(@(P, W) P .* repmat(W, size(P, 1), 1), ...
    conv_nb(s, x, y, z), repmat({comp_ws(s, x, y, z, m, b0)'}, nruns, 1), 'uni', false);
calc_br = @(s, x, y, z, m, b0) cellfun(@(S, M) (sum(S, 2) - M) .^ 2, ...
    pred_bs(s, x, y, z, m, b0), comp_bs(m, b0), 'uni', false);
calc_me = @(s, x, y, z, m, b0) sum(cell2mat(calc_br(s, x, y, z, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x(3), run_avgs, baseline);

end
