function flow = flow_out(v, t, in_flow)
% Calculates the output flow as a function of volume.
% 
% Adapted from code by JG (gru.stanford.edu/svn/matlab/balloonmodel.m)
% AS 9/2017

global which_tau tau_n tau_p tau alpha;

vt = dvdt(t, v, in_flow);
if (which_tau == 1)
    flow = v ^ (1 / alpha) + tau_p * vt;
    tau = tau_p;
else
    flow = v ^ (1 / alpha) + tau_n * vt;
    tau = tau_n;
end

end