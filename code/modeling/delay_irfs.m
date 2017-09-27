function irf = delay_irfs(tau_d, fs)
% Derive the delay channel impulse response function.
% 
% INPUT
%   1) tau_d: time constant of delay IRF
%   2) fs: sampling rate of IRF (Hz)
% 
% OUTPUT
%   irf: delay IRF (sampled at fs Hz)
% 
% AS 9/2017

% filter parameters from Watson
k = 1.33;       % ratio of time constants for different response stages
t2 = k * tau_d; % time constant of inhibitory mechanism
n1 = 9;         % number of stages in exctatory mechanism
n2 = 10;        % number of stages in inhibitory mechanism

% generate filters
time = 0:4000; % time in ms
for t = time
    % excitatory filter
    f1(t + 1) = ((tau_d * factorial(n1 - 1)) ^ -1) * ((t / tau_d) ^ (n1 - 1)) * exp(-t / tau_d);
    % inhibitory filter
    f2(t + 1) = ((t2 * factorial(n2 - 1)) ^ -1) * ((t / t2) ^ (n2 - 1)) * exp(-t / t2);
end

% derive delay IRF
irfS = f1;
irf = f1 - f2;
% normalize max of S and T IRFs
irf = irf * (max(irfS) / max(irf));
irf = resample(irf, 1, 1000 / fs);

end
