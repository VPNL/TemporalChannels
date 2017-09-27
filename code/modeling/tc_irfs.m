function irf = tc_irfs(channel, tau, fs)
% Derive the sustained and transient impulse response functions based on
% the formulation proposed by Watson (1986).
% 
% INPUT
%   1) channel: 'S' (sustained), 'T' (transient), or 'D' (delay)
%   2) tau: time constant of gamma function (ms)
%   3) fs: sampling rate of IRF (Hz)
% 
% OUTPUT
%   irf: either the sustained or transient IRF (sampled at fs Hz)
% 
% AS 2/2017

channel = upper(channel(1));
if nargin < 2; tau = 4.93; end
if nargin < 3; fs = 1000; end

% filter parameters from Watson
k = 1.33;     % ratio of time constants for different response stages
t2 = k * tau; % time constant of inhibitory mechanism
n1 = 9;       % number of stages in exctatory mechanism
n2 = 10;      % number of stages in inhibitory mechanism

% generate filters
time = 0:1499; % time in ms
for t = time
    % excitatory filter
    f1(t + 1) = ((tau * factorial(n1 - 1)) ^ -1) * ((t / tau) ^ (n1 - 1)) * exp(-t / tau);
    % inhibitory filter
    f2(t + 1) = ((t2 * factorial(n2 - 1)) ^ -1) * ((t / t2) ^ (n2 - 1)) * exp(-t / t2);
end

% derive sustained and transient IRFs
irfS = f1;
irfT = f1 - f2;
% normalize max of S and T IRFs
irfT = irfT * (max(irfS) / max(irfT));

% output sustained or transient IRF
if strcmp(channel,'S')
    irf = resample(irfS, 1, 1000 / fs)';
elseif sum(strcmp(channel,{'T' 'D'}))
    irf = resample(irfT, 1, 1000 / fs)';
else
    error('Unexpected channel input.');
end

% clip values from end that are close to zero for efficiency
zclip = length(irf); zthresh = max(irf) / 1000;
while abs(irf(zclip)) < zthresh
    zclip = zclip - 1;
end
irf = irf(1:zclip);

end

