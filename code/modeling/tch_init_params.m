function [params, irfs] = tch_init_params(model_type, nsess, fs)
% Initializes model parameters for given type of model
% 
% INPUTS
%   1) model_type: descriptor for type of model initialize
%     'glm'     -- irfs = hrf
%     'htd'     -- irfs = hrf, dhrf
%     'balloon' -- params = tauMTT, tau_n, tau_p, E0, V0, alpha
%     'cts-pow' -- params = tau1, epsilon; irfs = nrfS, hrf
%     'cts-div' -- params = tau1, sigma; irfs = nrfS, hrf
%     'dcts'    -- params = tau1, tau2, sigma; irfs = lpf, nrfS, hrf
%     '2ch-lin-lin'  -- irfs = nrfS, nrfT, hrf
%     '2ch-lin-quad' -- irfs = nrfS, nrfT, hrf
%     '2ch-lin-rect' -- irfs = nrfS, nrfT, hrf
%     '2ch-pow-quad' -- params = epsilon; irfs = nrfS, nrfT, hrf
%     '2ch-pow-rect' -- params = epsilon; irfs = nrfS, nrfT, hrf
%     '2ch-div-quad' -- params = sigma; irfs = lpf, nrfS, nrfT, hrf
%     '2ch-opt'      -- params = tau_s, tau_t; irfs = nrfS, nrfT, hrf
%     '3ch-lin-quad' -- params = tau_d, kappa; irfs = nrfS, nrfT, nrfD, hrf 
%     '3ch-lin-rect' -- params = tau_d, kappa; irfs = nrfS, nrfT, nrfD, hrf 
%     '3ch-pow-quad' -- params = epsilon, tau_d, kappa; irfs = nrfS, nrfT, nrfD, hrf 
%     '3ch-pow-rect' -- params = epsilon, tau_d, kappa; irfs = nrfS, nrfT, nrfD, hrf 
%     '3ch-opt'      -- params = tau_s, tau_t, tau_d; irfs = nrfS, nrfT, nrfD, hrf
%   2) nsess: number of sessions to setup parameters for 
%   3) fs: sampling rate for impulse response functions (Hz)
% 
% OUTPUTS
%   1) params: structures of model parameters for each session
%   2) irfs: structure of model impulse response functions for each session 
% 
% AS 2/2017

if nargin ~= 3; error('Unexpected input arguements.'); end
hrf = canonical_hrf(1 / fs, [5 14 28]);
dhrf = [diff(hrf); 0]; dhrf = dhrf * (max(hrf) / max(dhrf));
% default paramters for CTS family of moddels
epsilon = 0.1; tau1 = 100; tau2 = 150; sigma = 0.1;
% default paramters for multi-channel moddels
tau_s = 4.93; tau_t = 4.93; tau_d = 4.93; kappa = 1.33;
% exponential time constants for decay
tau_de = 500; tau_ae = 30000;
% default paramters for balloon model (see Chen & Glover, 2015)
tau_p = 25; tau_n = 25; tauMTT = 2.5; alpha = 0.4; E0 = 0.4; V0 = 0.03; 
% parameters of gamma IRF for balloon model
tau_g = 0.24; lag = 1; exponent = 2;
% setup structs
params = struct; irfs = struct;

switch model_type
    case 'glm'
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '1ch-lin'
        irfs.hrf = repmat({hrf}, 1, nsess);
    case 'htd'
        irfs.hrf = repmat({hrf}, 1, nsess);
        irfs.dhrf = repmat({dhrf}, 1, nsess);
    case 'balloon'
        params.tau_p = tau_p;   % viscoelastic time constant for inflation
        params.tau_n = tau_n;   % viscoelastic time constant for deflation
        params.E0 = E0;         % resting oxygen extraction fraction
        params.V0 = V0;         % resting blood volume
        params.tauMTT = tauMTT; % transit time through balloon
        params.alpha = alpha;   % steady-state flow-volume relationship
        % time constants for gamma function
        params.delta_t = 1 / fs; t_lag = (0:1 / fs:40) - lag;
        gamma_n = ((t_lag / tau_g) .^ (exponent - 1) .* exp(-t_lag / tau_g));
        gamma_d = (tau_g * factorial(exponent - 1));
        irfs.gamma = repmat({rectify(gamma_n ./ gamma_d)}, 1, nsess);
        % constants for calculating signal at 3T
        params.k1 = 6.7; params.k2 = 2.73; params.k3 = 0.57;
    case '1ch-balloon'
        params.tau_p = tau_p;   % viscoelastic time constant for inflation
        params.tau_n = tau_n;   % viscoelastic time constant for deflation
        params.E0 = E0;         % resting oxygen extraction fraction
        params.V0 = V0;         % resting blood volume
        params.tauMTT = tauMTT; % transit time through balloon
        params.alpha = alpha;   % steady-state flow-volume relationship
        % time constants for gamma function
        params.delta_t = 1 / fs; t_lag = (0:1 / fs:40) - lag;
        gamma_n = ((t_lag / tau_g) .^ (exponent - 1) .* exp(-t_lag / tau_g));
        gamma_d = (tau_g * factorial(exponent - 1));
        irfs.gamma = repmat({rectify(gamma_n ./ gamma_d)}, 1, nsess);
        % constants for calculating signal at 3T
        params.k1 = 6.7; params.k2 = 2.73; params.k3 = 0.57;
    case 'cts-pow'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau1 = repmat({tau1}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '1ch-pow'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau1 = repmat({tau1}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case 'cts-div'
        params.tau1 = repmat({tau1}, 1, nsess);
        params.sigma = repmat({sigma}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '1ch-div'
        params.tau1 = repmat({tau1}, 1, nsess);
        params.sigma = repmat({sigma}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case 'dcts'
        params.tau1 = repmat({tau1}, 1, nsess);
        params.tau2 = repmat({tau2}, 1, nsess);
        params.sigma = repmat({sigma}, 1, nsess);
        lpf = exp(-(0:999) / tau2); lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '1ch-dcts'
        params.tau1 = repmat({tau1}, 1, nsess);
        params.tau2 = repmat({tau2}, 1, nsess);
        params.sigma = repmat({sigma}, 1, nsess);
        lpf = exp(-(0:999) / tau2); lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '1ch-exp'
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
    case '1ch-cexp'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
    case '2ch-lin-htd'
        irfs.hrf = repmat({hrf}, 1, nsess);
        irfs.dhrf = repmat({dhrf}, 1, nsess);
    case '2ch-lin-lin'
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-lin-quad'
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-lin-rect'
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-pow-quad'
        params.epsilon = repmat({epsilon}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-pow-rect'
        params.epsilon = repmat({epsilon}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-div-quad'
        params.sigma = repmat({sigma}, 1, nsess);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-dcts-quad'
        params.tau2 = repmat({tau2}, 1, nsess);
        params.sigma = repmat({0.1}, 1, nsess);
        lpf = exp(-(0:999) / tau2); lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsess);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-opt'
        params.tau_s = repmat({tau_s}, 1, nsess);
        params.tau_t = repmat({tau_t}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-exp-quad'
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
    case '2ch-exp-rect'
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
    case '2ch-cexp-quad'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
    case '2ch-cexp-rect'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
    case '3ch-lin-quad'
        params.tau_d = repmat({tau_d}, 1, nsess);
        params.kappa = repmat({kappa}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        nrfD = tch_irfs('D', tau_d, kappa, fs);
        irfs.nrfD = repmat({nrfD}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '3ch-lin-rect'
        params.tau_d = repmat({tau_d}, 1, nsess);
        params.kappa = repmat({kappa}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        nrfD = tch_irfs('D', tau_d, kappa, fs);
        irfs.nrfD = repmat({nrfD}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '3ch-pow-quad'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_d = repmat({tau_d}, 1, nsess);
        params.kappa = repmat({kappa}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        nrfD = tch_irfs('D', tau_d, kappa, fs);
        irfs.nrfD = repmat({nrfD}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '3ch-pow-rect'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_d = repmat({tau_d}, 1, nsess);
        params.kappa = repmat({kappa}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        nrfD = tch_irfs('D', tau_d, kappa, fs);
        irfs.nrfD = repmat({nrfD}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '3ch-opt'
        params.tau_s = repmat({tau_s}, 1, nsess);
        params.tau_t = repmat({tau_t}, 1, nsess);
        params.tau_d = repmat({tau_d}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        nrfD = tch_irfs('D', tau_d, kappa, fs);
        irfs.nrfD = repmat({nrfD}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '3ch-lin-quad-exp'
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
    case '3ch-lin-rect-exp'
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
    case '3ch-pow-quad-exp'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
    case '3ch-exp-quad-exp'
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
    case '3ch-exp-rect-exp'
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
    case '3ch-cexp-quad-exp'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
    case '3ch-cexp-rect-exp'
        params.epsilon = repmat({epsilon}, 1, nsess);
        params.tau_ae = repmat({tau_ae}, 1, nsess);
        params.tau_de = repmat({tau_de}, 1, nsess);
        nrfS = tch_irfs('S', tau_s, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = tch_irfs('T', tau_t, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsess);
        delay_exp = exp(-(1:12000) / tau_de);
        irfs.delay_exp = repmat({delay_exp}, 1, nsess);
end

end
