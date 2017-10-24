function obj_fun = tch_obj_fun_3ch_exp_quad_exp(roi, model)
% Generates anonymous objective function that can be passed to fmincon for
% the 3ch-exp-quad-exp model (3-channel model with adapted sustained,
% quadratic transient, and exponential delay channels).
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

if ~strcmp(model.type, '3ch-exp-quad-exp'); error('Incompatible model type'); end
stim = model.stim; nruns = size(stim, 1); npreds = size(stim{1}, 2);
irfs = model.irfs; fs = model.fs; tr = model.tr;
run_avgs = roi.run_avgs; baseline = roi.baseline;
param_names = fieldnames(model.params); nparams = length(param_names);

adapt_fun = @(y) exp(-(1:60000) / (y * 1000));
delay_fun = @(z) exp(-(1:12000) / (z * 1000));
conv_snS = @(x, y) cellfun(@(X, Y, ON, OFF) code_exp_decay(X, ON, OFF, Y, fs), ...
    cellfun(@(XX) convolve_vecs(XX, irfs.nrfS{1}, 1, 1), x, 'uni', false), ...
    repmat({adapt_fun(y)}, nruns, 1), model.onsets, model.offsets, 'uni', false);
conv_snT = @(x) cellfun(@(X) convolve_vecs(X, irfs.nrfT{1}, 1, 1) .^ 2, ...
    x, 'uni', false);
doffsets = cellfun(@(X, Y) [X(2:end) Y], model.onsets, model.run_durs, 'uni', false);
conv_snD = @(x, z) cellfun(@(X, Z, ON, OFF) code_exp_decay(X, ON, OFF, Z, fs), ...
    cellfun(@code_delay_act, x, 'uni', false), repmat({delay_fun(z)}, nruns, 1), ...
    model.offsets, doffsets, 'uni', false);
conv_nbS = @(x, y) cellfun(@(NS) convolve_vecs(NS, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snS(x, y), 'uni', false);
conv_nbT = @(x) cellfun(@(NT) convolve_vecs(NT, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snT(x), 'uni', false);
conv_nbD = @(x, z) cellfun(@(ND) convolve_vecs(ND, irfs.hrf{1}, fs, 1 / tr), ...
    conv_snD(x, z), 'uni', false);
pred_bsS = @(x, y, b) cellfun(@(PS, BS) PS .* repmat(BS, size(PS, 1), 1), ...
    conv_nbS(x, y), repmat({b}, nruns, 1), 'uni', false);
pred_bsT = @(x, b) cellfun(@(PT, BT) PT .* repmat(BT, size(PT, 1), 1), ...
    conv_nbT(x), repmat({b}, nruns, 1), 'uni', false);
pred_bsD = @(x, z, b) cellfun(@(PD, BD) PD .* repmat(BD, size(PD, 1), 1), ...
    conv_nbD(x, z), repmat({b}, nruns, 1), 'uni', false);
comp_bs = @(m, b0) cellfun(@(M, B0) M - repmat(B0, size(M, 1), 1), ...
    m, b0, 'uni', false);
calc_br = @(x, y, z, b, m, b0) cellfun(@(SS, ST, SD, M) (sum([SS ST SD], 2) - M) .^ 2, ...
    pred_bsS(x, y, b(1:npreds)), pred_bsT(x, b([1:npreds] + npreds)), ...
    pred_bsD(x, z, b([1:npreds] + 2 * npreds)), comp_bs(m, b0), 'uni', false);
calc_me = @(x, y, z, b, m, b0) sum(cell2mat(calc_br(x, y, z, b, m, b0)));
obj_fun = @(x) calc_me(stim, x(1), x(2), x([1:npreds * 3] + nparams), run_avgs, baseline);

end
