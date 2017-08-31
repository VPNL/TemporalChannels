function blims = axesTS(name, xlab, ylab, lims)
% Draws custom 2ch square axes for plotting relative channel contributions.
% 
% Inputs
%   1) name: title for plot
%   2) xlab: x-axis label
%   3) ylab: y-axis label
%   4) lims: limits of plotted axes
% 
% Outputs
%   1) blims: limits of hidden axes
% 
% AS 8/2017

% find limits of plotting box
blims = lims; shift = range(lims) * 0.04; blims(1) = lims(1) - shift;
tl = range(blims) * 0.04; blims(1) = blims(1) - tl * 5.5;
xlim(blims); ylim(blims);

% draw title and axis labels
text(mean(lims), lims(2), name, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 7, 'FontName', 'Helvetica', 'FontWeight', 'bold');

text(mean(lims), blims(1), xlab, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 7, 'FontName', 'Helvetica');

text(blims(1), mean(lims), ylab, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 7, 'FontName', 'Helvetica', 'Rotation', 90);

% plot axis lines for origin and ruler
plot([0 lims(2)], [0 0], '--', 'Color', [.8 .8 .8], 'LineWidth', 0.5);
plot([0 lims(2)], [0 0] - shift, 'k', 'LineWidth', 0.5);
plot([0 0], [0 lims(2)], '--', 'Color', [.8 .8 .8], 'LineWidth', 0.5);
plot([0 0] - shift, [0 lims(2)], 'k', 'LineWidth', 0.5);

% plot major tick marks for both axes
tmajor = 0:lims(2);
for tt = tmajor
    plot([tt tt], [0 - shift 0 - shift - tl], 'k', 'LineWidth', 0.5);
    text(tt, lims(1) - shift - tl * 1.25, num2str(tt), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 6, 'FontName', 'Helvetica');
    plot([0 - shift 0 - shift - tl], [tt tt], 'k', 'LineWidth', 0.5);
    text(lims(1) - shift - tl * 1.25, tt, num2str(tt), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', 6, 'FontName', 'Helvetica');
end

% plot minor tick marks for both axes
tminor = 0:0.25:lims(2); tminor(rem(tminor, 1) == 0) = [];
for tt = tminor
    plot([tt tt], [0 - shift 0 - shift - tl / 2], 'k', 'LineWidth', 0.5);
    plot([0 - shift 0 - shift - tl / 2], [tt tt], 'k', 'LineWidth', 0.5);
end

% limit axis ranges and hide automatic axes
xlim(blims + tl); ylim(blims + tl); axis square; axis off;

end
