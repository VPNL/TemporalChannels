function [cx, cy] = draw_cross(x, y, w, h, col, lw)
% Plots an cross centered at (x, y) of width 2w and height 2h
% INPUTS
%   1) x: x-axis value of octogon origin
%   2) y: y-axis value of octogon origin
%   3) w: width from origin to horizontal extremes
%   4) h: height from origin to vertical extremes
%   5) col: octogon fill color
%   6) lw: line width
% 
% OUTPUT
%   1) cx: handle to octogon fill object
%   2) cy: handle to octogon fill object
% 
% AS 8/2017

if nargin < 5
    col = [0 0 0];
end
if nargin < 6
    lw = 3;
end

cx = plot([x - w x + w], [y y], 'Color', col, 'LineWidth', lw);
cy = plot([x x], [y - h y + h], 'Color', col, 'LineWidth', lw);

end
