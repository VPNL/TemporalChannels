function [params, irfs] = tch_init_params(model_type, nsessions, fs)
% Initializes model parameters for given type of model
% 
% INPUTS
%   1) model_type: descriptor for type of model initialize
%   2) nsessions: number of sessions to setup parameters for 
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
[tau_s, tau_t, tau_p] = deal(4.93); n1 = 9; n2 = 10; kappa = 1.33;
% exponential time constants for decay
tau_pe = 500; tau_ae = 20000;
% default paramters for balloon model (see Chen & Glover, 2015)
tauP = 25; tauN = 25; tauMTT = 2.5; alpha = 0.4; E0 = 0.4; V0 = 0.03; 
% parameters of gamma IRF for balloon model
tau_g = 0.24; lag = 1; exponent = 2;
% setup structs
params = struct; irfs = struct;

switch model_type
    case '1ch-lin'
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '1ch-balloon'
        params.tauP = tauP;   % viscoelastic time constant for inflation
        params.tauN = tauN;   % viscoelastic time constant for deflation
        params.E0 = E0;         % resting oxygen extraction fraction
        params.V0 = V0;         % resting blood volume
        params.tauMTT = tauMTT; % transit time through balloon
        params.alpha = alpha;   % steady-state flow-volume relationship
        % time constants for gamma function
        params.delta_t = 1 / fs; t_lag = (0:1 / fs:40) - lag;
        gamma_n = ((t_lag / tau_g) .^ (exponent - 1) .* exp(-t_lag / tau_g));
        gamma_d = (tau_g * factorial(exponent - 1));
        irfs.gamma = repmat({rectify(gamma_n ./ gamma_d)}, 1, nsessions);
        % constants for calculating signal at 3T
        params.k1 = 6.7; params.k2 = 2.73; params.k3 = 0.57;
    case '1ch-pow'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau1 = repmat({tau1}, 1, nsessions);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '1ch-div'
        params.tau1 = repmat({tau1}, 1, nsessions);
        params.sigma = repmat({sigma}, 1, nsessions);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '1ch-dcts'
        params.tau1 = repmat({tau1}, 1, nsessions);
        params.tau2 = repmat({tau2}, 1, nsessions);
        params.sigma = repmat({sigma}, 1, nsessions);
        lpf = exp(-(0:999) / tau2); lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsessions);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '1ch-exp'
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
    case '1ch-cexp'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
    case '1ch-quad-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-htd'
        irfs.hrf = repmat({hrf}, 1, nsessions);
        irfs.dhrf = repmat({dhrf}, 1, nsessions);
    case '2ch-lin-lin'
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-quad'
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        %nrfT = -tch_irfs('T', tau_t, 9, 17, .5);
        %nrfT = -tch_irfs('T', tau_t, 5, 9, .5);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-rect'
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-pow-quad'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-pow-rect'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-div-quad'
        params.sigma = repmat({sigma}, 1, nsessions);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-dcts-quad'
        params.tau2 = repmat({tau2}, 1, nsessions);
        params.sigma = repmat({0.1}, 1, nsessions);
        lpf = exp(-(0:999) / tau2); lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsessions);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-exp-quad'
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
    case '2ch-exp-rect'
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
    case '2ch-cexp-quad'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
    case '2ch-cexp-rect'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
    case '3ch-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_t = repmat({tau_t}, 1, nsessions);
        params.tau_d = repmat({tau_p}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        nrfP = tch_irfs('P', tau_p, n1, n2, kappa, fs);
        irfs.nrfP = repmat({nrfP}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '3ch-lin-quad-exp'
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-lin-rect-exp'
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-pow-quad-exp'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-pow-rect-exp'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-exp-quad-exp'
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-exp-rect-exp'
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-cexp-quad-exp'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '3ch-cexp-rect-exp'
        params.epsilon = repmat({epsilon}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
    case '2ch-lin-quad-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.n1 = repmat({n1}, 1, nsessions);
        params.n2 = repmat({n2}, 1, nsessions);
        params.kappa = repmat({kappa}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-rect-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.n1 = repmat({n1}, 1, nsessions);
        params.n2 = repmat({n2}, 1, nsessions);
        params.kappa = repmat({kappa}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '3ch-lin-quad-exp-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.n1 = repmat({n1}, 1, nsessions);
        params.n2 = repmat({n2}, 1, nsessions);
        params.kappa = repmat({kappa}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '3ch-lin-rect-exp-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.n1 = repmat({n1}, 1, nsessions);
        params.n2 = repmat({n2}, 1, nsessions);
        params.kappa = repmat({kappa}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '1ch-exp-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-exp-quad-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '3ch-exp-quad-exp-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau_pe = repmat({tau_pe}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        persist_exp = exp(-(1:12000) / tau_pe);
        irfs.persist_exp = repmat({persist_exp}, 1, nsessions);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-exp-cquad-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau2 = repmat({tau2}, 1, nsessions);
        params.sigma = repmat({sigma}, 1, nsessions);
        lpf = exp(-(0:999) / tau2); lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '3ch-exp-quad-crect-opt'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.tau_p = repmat({tau_p}, 1, nsessions);
        params.epsilon = repmat({epsilon}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        nrfP = tch_irfs('P', tau_p, n1, n2, kappa, fs);
        irfs.nrfP = repmat({nrfP}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-exp-crect'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.epsilon = repmat({epsilon}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-crect'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.epsilon = repmat({epsilon}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-exp-cquad'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_t = repmat({tau_t}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_t, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-cquad'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.epsilon = repmat({epsilon}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-exp-dquad'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.tau_ae = repmat({tau_ae}, 1, nsessions);
        params.sigma = repmat({sigma}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        adapt_exp = exp(-(1:60000) / tau_ae);
        irfs.adapt_exp = repmat({adapt_exp}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
    case '2ch-lin-dquad'
        params.tau_s = repmat({tau_s}, 1, nsessions);
        params.sigma = repmat({sigma}, 1, nsessions);
        nrfS = tch_irfs('S', tau_s, n1, n2, kappa, fs);
        irfs.nrfS = repmat({nrfS}, 1, nsessions);
        nrfT = tch_irfs('T', tau_s, n1, n2, kappa, fs);
        irfs.nrfT = repmat({nrfT}, 1, nsessions);
        irfs.hrf = repmat({hrf}, 1, nsessions);
end

end
