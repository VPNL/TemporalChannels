function tch_axes(names, xlab, ylab, xlims, ylims)
% Draws custom axes for plotting graphs with shared y-axis.
% 
% Inputs
%   1) names: title for each subplot
%   2) xlab: x-axis label
%   3) ylab: y-axis label
%   4) xlims: limits of plotted x-axes (N by 2 matrix)
%   5) ylims: limits of plotted y-axis (1 x 2 vector)
% 
% AS 8/2017

% find overall range of plotting area and format axes
names = force_cell(names); xrange = [min([min(xlims(:)) 0]) max(xlims(:))];
xlim(xrange); ylim([ylims(1) ylims(2)]); hold on;
set(gca, 'XColor', 'w', 'TickDir', 'out', 'YMinorTick', 'on', 'FontSize', 7);

% set title, x-axes, and axis labels
for xx = 1:size(xlims, 1)
    text(mean(xlims(xx, :)), ylims(2), names{xx}, ...
        'FontSize', 7, 'FontWeight', 'bold', 'FontName', 'Helvetica', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top');
    plot(xlims(xx, :), [0 0], 'k-', 'LineWidth', 0.5);
end
xlab_x = mean([min(xlims(:)) max(xlims(:))]);
xlab_y = get(gca, 'XLabel'); xlab_y = xlab_y.Position(2);
text(xlab_x, xlab_y, xlab, 'FontSize', 7, 'FontName', 'Helvetica', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
ylabel(ylab, 'FontSize', 7);

end
