function tc_axes(name, xlab, ylab, xlims, ylims)
% Draws custom axes for plotting peristimulus time series.
% 
% Inputs
%   1) name: title for plot
%   2) xlab: x-axis label
%   3) ylab: y-axis label
%   4) xlims: limits of plotted x-axes (vector or matrix)
%   5) ylims: limits of plotted y-axis (vector)
% 
% AS 8/2017

% find overall range of plotting area and format axes
xrange = [min([min(xlims(:)) 0]) max(xlims(:))];
xlim(xrange); ylim([ylims(1) ylims(2)]); hold on;
set(gca, 'XColor', 'w', 'FontSize', 6, 'FontName', 'Helvetica', ...
    'TickDir', 'out', 'YMinorTick', 'on');

% set title, x-axes, and axis labels
title(name, 'FontSize', 7, 'FontName', 'Helvetica');
ylabel(ylab, 'FontSize', 7, 'FontName', 'Helvetica');
for xx = 1:size(xlims)
    plot(xlims(xx, :), [0 0], 'k-', 'LineWidth', 0.5);
end
xlab_x = mean([min(xlims(:)) max(xlims(:))]);
xlab_y = get(gca, 'XLabel'); xlab_y = xlab_y.Position(2);
text(xlab_x, xlab_y, xlab, 'FontSize', 7, 'FontName', 'Helvetica', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');

end
