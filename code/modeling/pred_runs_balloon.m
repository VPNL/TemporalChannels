function model = pred_runs_balloon(t, stim, varargin)
% Generates run predictors using the standard linear model.

% setup some constants
tau_p = 20; tau_n = 20; lag = 1; tau = 0.24;
ex = 2; delta_t = diff(t); delta_t = delta_t(1);
E0 = 0.4;    % resting oxygen extraction rate
V0 = 0.04;   % venous blood volume fraction
tau_m = 3;   % time of the volume/deoxyhemoglobin variables
wt = 1;      % begin with positive tau
alpha = 0.4;

% store values in params struct
params.tau_p = tau_p; params.tau_n = tau_n; params.tau = tau;
params.lag = lag; params.delta_t = delta_t; params.E0 = E0; params.V0 = V0;
params.tau_m = tau_m; params.wt = wt; params.alpha = alpha;

% convolve input with a gamma function
gamma_n = ((t - lag) / tau) .^ (ex - 1) .* exp(-(t - lag) / tau);
gamma_d = tau * factorial(ex - 1);
gamma_fun = gamma_n ./ gamma_d;
in_flow = conv(model.stim, gamma_fun) / sum(gamma_fun);
in_flow = reshape(.7 * in_flow + 1, 1, []);

% simulate model
model = simulateBalloonModel(t, in_flow, params);

% set initial values
[v, q, IN_FLOW, OUT_FLOW, CMRO2] = deal(0); S = 0; OEF = E0;

% get the simulated values of all variables
for i = 1:length(t)
    v(i + 1) = runge_kutta(delta_t, @dvdt, t(i), v(i), in_flow(i), params);
    OUT_FLOW(i + 1) = out_flow(v(i + 1), t(i), in_flow(i), params);
    q(i + 1) = runge_kutta(delta_t, @dqdt, t(i), q(i), v(i), in_flow(i), params);
    S(i + 1) = sig(q(i + 1), v(i + 1), E0, V0);
    OEF(i + 1) = E(in_flow(i), E0);
    CMRO2(i + 1) = (OEF(i + 1)/E0) * IN_FLOW(i);
    IN_FLOW(i + 1) = in_flow(i);
end
% set time variable to include initial time point
t = [0 t];

% pack it up
model.t = t;
model.S = S;
model.IN_FLOW = IN_FLOW;
model.OUT_FLOW = OUT_FLOW;
model.q = q;
model.v = v;
model.OEF = OEF;
model.CMRO2 = CMRO2;

end


%% Rate of change in normalized deoxyhemoglobin (Equation 5.2)
function out_val = dqdt(t, q, v, in_flow, params)

tau_m = params.tau_m;
E0 = params.E0;

if (iscell(v))
    in_flow = v{2};
    v = v{1};
end

POS_FLOW = in_flow * (E(in_flow, E0) / E0);
NEG_FLOW = out_flow(v, t, in_flow, params) * (q / v);
out_val = (1 / tau_m) * (POS_FLOW - NEG_FLOW);

end

%% rate of change in volume as a function of volume (Equation 5.2)
function [out_val, wt] = dvdt(t, v, in_flow, params)

tau_m = params.tau_m;
wt = params.wt;
alpha = params.alpha;
tau_p = params.tau_p;
tau_n = params.tau_n;

% get fin if it is passed by Runge-Kutta
if (iscell(in_flow))
    in_flow = in_flow{1};
end

if (wt == 1)
    POS_FLOW = ((1 / tau_m) * (in_flow - v ^ (1 / alpha)));
    out_val = POS_FLOW / (1 + (tau_p / tau_m));
else
    NEG_FLOW = ((1 / tau_m) * (in_flow - v ^ (1 / alpha)));
    out_val = NEG_FLOW / (1 + (tau_n / tau_m));
end

if out_val < 0
    wt = 2;
else
    wt = 1;
end

end


%% Function that computes the output flow as a function of volume 
function [out_val, tau] = out_flow(v, t, in_flow, params)

wt = params.wt;
alpha = params.alpha;
tau_p = params.tau_p;
tau_n = params.tau_n;

if (wt == 1)
    out_val = v ^ (1 / alpha) + tau_p * dvdt(t, v, in_flow, params);
    tau = tau_p;
else
    out_val = v ^ (1 / alpha) + tau_n * dvdt(t, v, in_flow, params);
    tau = tau_n;
end

end

%% Function that determines signal strength as a function of:
% q -- the normalized total deoxyhemoglobin
% v -- the normalized volume
% The equation uses V0 the fraction of blood volume
% in veins to compute the intravascular and
% extravascular components of signal change. The
% constants are taken from Ogawa et al.'s monte
% carlo simulation and other studies and are based
% on a 1.5T B0 with TE = 40 ms.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out_val = sig(q, v, E0, V0)

k1 = 7 * E0;
k2 = 2;
k3 = 2 * E0 - 0.2;
out_val = V0 * (k1 * (1 - q) + k2 * (1 - (q / v)) + k3 * (1 - v));

end

%% Function that relates the Oxygen extraction rate to the flow. Eq (6)
function out_val = E(f, E0)

out_val = 1 - (1 - E0) .^ (1 ./ f);

end