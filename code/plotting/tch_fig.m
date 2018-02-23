function fig = tch_fig(name, pos)

if nargin < 1
    name = '';
end
if nargin < 2
    pos = [.1 .1 .5 .5];
else
    rectify(pos);
end
fig = figure('Name', name, 'Color', 'w', 'Units', 'norm', 'Position', pos);
tch_set_axes; hold on;

end