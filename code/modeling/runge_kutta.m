function out_val = runge_kutta(delta_t, flow_fun, t, y, varargin)
% Performs fourth-order Runge-Kutta method for balloon model optimization.
% AS 9/2017 adapted from code written by Justin Gardner

% k1 is the Euler step
if (nargin > 4)
    k1 = delta_t * feval(flow_fun, t, y, varargin);
else
    k1 = delta_t * feval(flow_fun, t, y);
end

% k2 is the step using the derivative at the midpoint of the euler step
if (nargin > 4)
    k2 = delta_t * feval(flow_fun, t + delta_t / 2, y + k1 / 2, varargin);
else
    k2 = delta_t * feval(flow_fun, t + delta_t / 2, y + k1 / 2);
end

% k3 is the step using the derivative at the midpoint of the step above
if (nargin > 4)
    k3 = delta_t * feval(flow_fun, t + delta_t / 2, y + k2 / 2, varargin);
else
    k3 = delta_t * feval(flow_fun, t + delta_t / 2, y + k2 / 2);
end

% k4 is the step using the derivative at the endpoint of the step above
if (nargin > 4)
    k4 = delta_t * feval(flow_fun, t + delta_t, y + k3, varargin);
else
    k4 = delta_t * feval(flow_fun, t + delta_t, y + k3);
end

% compute final output
out_val = y + k1 / 6 + k2 / 3 + k3 / 3 + k4 / 6;

end