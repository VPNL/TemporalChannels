function bl = axesTS(name, xlab, ylab, lims)
% Draws custom 2ch square axes for plotting relative channel contributions.
% AS 8/2017

% find limits of plotting box
bl = lims; shift = range(lims) * 0.02; bl(1) = lims(1) - shift;
tl = range(bl) * 0.02; bl(1) = bl(1) - tl * 6;
xlim(bl); ylim(bl);

% draw title and axis labels
text(mean(lims), lims(2), name, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 7, 'FontName', 'Helvetica', 'FontWeight', 'bold');

text(mean(lims), bl(1), xlab, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'FontSize', 7, 'FontName', 'Helvetica');

text(bl(1), mean(lims), ylab, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 7, 'FontName', 'Helvetica', 'Rotation', 90);

% plot axis lines for origin and ruler
plot([0 lims(2)], [0 0], '--', 'Color', [.8 .8 .8], 'LineWidth', 0.5);
plot([0 lims(2)], [0 0] - shift, 'k', 'LineWidth', 0.5);
plot([0 0], [0 lims(2)], '--', 'Color', [.8 .8 .8], 'LineWidth', 0.5);
plot([0 0] - shift, [0 lims(2)], 'k', 'LineWidth', 0.5);

% plot tick marks for both axes
for tt = 0:lims(2)
    plot([tt tt], [0 - shift 0 - shift - tl], 'k', 'LineWidth', 0.5);
    text(tt, lims(1) - shift - tl * 1.25, num2str(tt), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 6, 'FontName', 'Helvetica');
    plot([0 - shift 0 - shift - tl], [tt tt], 'k', 'LineWidth', 0.5);
    text(lims(1) - shift - tl * 1.25, tt, num2str(tt), ...
        'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
        'FontSize', 6, 'FontName', 'Helvetica');
end

% limit axis ranges and hide automatic axes
xlim(bl + tl); ylim(bl + tl); axis square; axis off;

end
