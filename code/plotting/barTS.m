function [ymin, ymax] = barTS(y, col, x1)
% Generates clean version of simple bar graph with error bars.
% 
% INPUTS
%   1) y: matrix of data (N x P)
%   2) col: plotting color for bar
%   3) x1: x-axis location of first bar in series
% 
% OUTPUTS
%   1) ymin: minimum value of bar - err
%   2) ymax: maximum value of bar + err
% 
% AS 5/2017

% check inputs
if nargin < 2 || isempty(col); col = [0 0 0]; end
if size(col, 1) == 1; col = repmat(col, size(y, 2), 1); end
if nargin < 3; x1 = 0; end

% compute mean and standard error of y
y_m = mean(y, 1);
if size(y, 1) > 1
    err = std(y) / sqrt(size(y, 1) - 1);
else
    err = zeros(size(y));
end

% plot bars and error bars
xv = [-1 1 1 -1]; xve = 0.1 * xv + x1;  yv = [0 0 1 1]; yve = [-1 -1 1 1];
xvec = 0.45 * xv + x1; xcnt = 0; 
for bb = 1:length(y_m)
    yvec = y_m(bb) * yv; xcnt = xcnt + 1;
    patch(xcnt + xvec, yvec, col(bb, :));
    yvec_e = repmat(y_m(bb), 1, 4) + yve * err(bb);
    patch(xcnt + xve, yvec_e, [1 1 1]);
end

% get limits of y-axis
ymin = min([min(y_m - err) 0]); ymax = max(y_m + err);

end
