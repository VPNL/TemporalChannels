function [params, irfs] = init_params(model_type, nsess, fs)
% Initializes model parameters for given type of model
% 
% INPUTS
%   1) model_type: descriptor for type of model initialize
%     'standard' -- irfs = hrf
%     'htd'      -- irfs = hrf, dhrf
%     '2ch-lin'  -- irfs = nrfS, nrfT, hrf
%     'cts-pow'  -- params = tau1, epsilon; irfs = nrfS, hrf
%     'cts-div'  -- params = tau1, sigma; irfs = nrfS, hrf
%     'dcts'     -- params = tau1, tau2, sigma; irfs = lpf, nrfS, hrf
%     '2ch'      -- irfs = nrfS, nrfT, hrf
%     '2ch-rect' -- irfs = nrfS, nrfT, hrf
%     '2ch-pow'  -- params = epsilon; irfs = nrfS, nrfT, hrf
%     '2ch-div'  -- params = sigma; irfs = lpf, nrfS, nrfT, hrf
%     '2ch-dcts' -- params = tau2, sigma; irfs = lpf, nrfS, nrfT, hrf
%     '2ch-opt' -- params = tau1, tau2, sigma; irfs = lpf, nrfS, nrfT, hrf
%   2) nsess: number of sessions to setup parameters for 
%   3) fs: sampling rate for impulse response functions (Hz)
% 
% OUTPUTS
%   1) params: structures of model parameters for each session
%   2) irfs: structure of model impulse response functions for each session 
% 
% AS 2/2017

if nargin ~= 3
    error('Unexpected input arguements.');
end
hrf = canonical_hrf(1 / fs, [5 14 28]);
dhrf = [diff(hrf); 0]; dhrf = dhrf * (max(hrf) / max(dhrf));
epsilon = 0.1; tau1 = 100; tau2 = 150; sigma = 0.1;
params = struct; irfs = struct;

switch model_type
    case 'standard'
        irfs.hrf = repmat({hrf}, 1, nsess);
    case 'htd'
        irfs.hrf = repmat({hrf}, 1, nsess);
        irfs.dhrf = repmat({dhrf}, 1, nsess);
    case '2ch-lin'
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case 'cts-pow'
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
    case 'dcts'
        params.tau1 = repmat({tau1}, 1, nsess);
        params.tau2 = repmat({tau2}, 1, nsess);
        params.sigma = repmat({sigma}, 1, nsess);
        lpf = exp(-(0:999) / tau2);
        lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch'
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-rect'
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-pow'
        params.epsilon = repmat({epsilon}, 1, nsess);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-div'
        params.sigma = repmat({sigma}, 1, nsess);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-dcts'
        params.tau2 = repmat({tau2}, 1, nsess);
        params.sigma = repmat({0.1}, 1, nsess);
        lpf = exp(-(0:999) / tau2);
        lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsess);
        nrfS = watson_irfs('S', fs);
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
    case '2ch-opt'
        params.tau1 = repmat({tau1}, 1, nsess);
        params.tau2 = repmat({tau2}, 1, nsess);
        params.sigma = repmat({0.1}, 1, nsess);
        lpf = exp(-(0:999) / tau2);
        lpf = lpf / sum(lpf);
        irfs.lpf = repmat({lpf}, 1, nsess);
        nrfS = (0:999) .* exp(-(0:999) / tau1);
        nrfS = resample(nrfS/sum(nrfS), 1, 1000 / fs)';
        irfs.nrfS = repmat({nrfS}, 1, nsess);
        nrfT = watson_irfs('T', fs);
        irfs.nrfT = repmat({nrfT}, 1, nsess);
        irfs.hrf = repmat({hrf}, 1, nsess);
end

end
