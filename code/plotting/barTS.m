function [ymin, ymax] = barTS(y, col)
% Generates clean version of simple bar graph with error bars.
% 
% INPUTS
%   1) y: matrix of data (N x P)
%   2) col: plotting color for bar
% 
% OUTPUTS
%   1) ymin: minimum value of bar - err
%   2) ymax: maximum value of bar + err
% 
% AS 5/2017

% set plotting color for each bar
if size(col, 1) == 1
    col = repmat(col, size(y, 2), 1);
end

% compute mean and standard error of y
y_m = mean(y, 1);
if size(y, 1) > 1
    err = std(y) / sqrt(size(y, 1) - 1);
else
    err = zeros(size(y));
end

% plot bars and error bars
vv_x = [-1 1 1 -1]; xvec = .45*vv_x; xvec_e = .1*vv_x;
vv_y = [0 0 1 1]; vv_ye = [-1 -1 1 1]; xcnt = 1; 
for bb = 1:length(y_m)
    yvec = y_m(bb) * vv_y;
    patch(xcnt + xvec, yvec, col(bb, :));
    yvec_e = repmat(y_m(bb), 1, 4) + vv_ye * err(bb);
    patch(xcnt + xvec_e, yvec_e, [1 1 1]);
    xcnt = xcnt+1;
end

% format
xlim([0 xcnt]);
ymin = min([floor(min(y_m - err)) 0]);
ymax = ceil(max(y_m + err));
ylim([ymin ymax]);
set(gca, 'FontSize', 8, 'XTick', 1:xcnt - 1, 'TickDir', 'out');

end
