function irf = delay_irf(fs)
% Derive biphasic impulse response function for delay activity
% 
% INPUT
%   fs: sampling rate of IRF in Hz
% 
% OUTPUT
%   irf: delay activity IRF (sampled at fs Hz)
% 
% AS 4/2017

% filter parameters 
k = 1.33;    % ratio of time constants for different response stages
t1 = 30;     % time constant of excitatory mechanism
t2 = k * t1; % time constant of inhibitory mechanism
n1 = 9;      % number of stages in exctatory mechanism
n2 = 10;     % number of stages in inhibitory mechanism

% generate filters
time = 0:999; % time in ms
for t = time
    % excitatory filter
    f1(t+1) = ((t1 * factorial(n1 - 1)) ^ -1) * ((t / t1) ^ (n1 - 1)) * exp(-t / t1);
    % inhibitory filter
    f2(t+1) = ((t2 * factorial(n2 - 1)) ^ -1) * ((t / t2) ^ (n2 - 1)) * exp(-t / t2);
end

% derive the IRF
irf = resample(f1 - f2, 1, 1000 / fs);

end
