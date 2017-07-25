function oct = draw_rect(x, y, w, h, fill_col, line_col)
% Plots a rectangle centered at (x, y) of width 2w and height 2h
% INPUTS
%   1) x: x-axis value of rectangle origin
%   2) y: y-axis value of rectangle origin
%   3) w: width from origin to horizontal extremes
%   4) h: height from origin to vertical extremes
%   5) fill_col: rectangle fill color
%   6) line_col: outline color
% 
% OUTPUT
%   1) oct: handle to octogon fill object
% 
% AS 6/2017

xx = [x - w; x + w; x + w; x - w]';
yy = [y - h; y - h; y + h; y + h]';
ff = fill(xx, yy, fill_col, 'EdgeColor', line_col);

end
