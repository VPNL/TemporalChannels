function [ymin, ymax] = tch_plot_ceil(x, y)
% Plots error bars for showing confidence intervals and noise ceilings.
% 
% INPUTS
%   1) x: range of x-axis to plot bar
%   2) y: matrix of data (N x P)
% 
% OUTPUTS
%   1) ymin: minimum value of bar - err
%   2) ymax: maximum value of bar + err
% 
% AS 7/2017

% compute mean and standard error of y
y_m = mean(y, 1);
if size(y, 1) > 1
    err = std(y) / sqrt(size(y, 1) - 1);
else
    err = zeros(size(y));
end

% plot mean and error region
vv_y = [0 0 1 1]; vv_ye = [-1 -1 1 1];
yvec_e = repmat(y_m, 1, 4) + vv_ye * err;
patch([x fliplr(x)], yvec_e, [.9 .9 .9], 'EdgeColor', [.9 .9 .9]);
plot(x, [y_m y_m], 'k:');

% format
ymin = min([floor(min(y_m - err)) 0]);
ymax = ceil(max(y_m + err));
ylim([ymin ymax]);

end
