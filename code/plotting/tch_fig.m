function fig = tch_fig(name, pos)

fig = figure('Name', name, 'Color', 'w', 'Units', 'norm', 'Position', pos);
tch_set_axes(gca); hold on;

end