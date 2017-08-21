function f = plot_betas(roi)
% Plots transient vs. sustained contributions for each stimulus category. 
% AS 8/2017

roi = roi(1);
betas = cell2mat(roi.model.betas');
means = mean(betas);
sems = std(betas) / sqrt(size(betas, 1) - 1);
ncats = length(roi.model.cat_list);
fill_cols = [.5 .5 1; 1 .5 .5; .5 .5 .5];
line_cols = [0 0 1; 1 0 0; 0 0 0];
xm = means(:, ncats + 1:ncats * 2);
ym = means(:, 1:ncats);
xe = sems(:, ncats + 1:ncats * 2);
ye = sems(:, 1:ncats);

xmax = max(xm + xe); ymax = max(ym + ye);
xmin = min(xm - xe); ymin = min(ym - ye);
lims = [min([0 xmin ymin]) max([ceil(xmax) ceil(ymax)])];

f = figTS([roi.model.type ' contributions'], [.1 .1 .3 .5]);
bl = axesTS(roi.nickname, 'Transient \beta_T', 'Sustained \beta_S', lims);
br = range(bl); lcnt = bl(2) - br * 0.05; xl = br * 0.05;
for cc = 1:ncats
    oo = draw_oct(xm(cc), ym(cc), xe(cc), ye(cc), fill_cols(cc, :), line_cols(cc, :));
    oo = draw_oct(xl, lcnt, br * 0.02, br * 0.02, fill_cols(cc, :), line_cols(cc, :));
    text(xl + br * 0.03, lcnt, roi.model.cat_list{cc}, ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', 6, 'FontName', 'Helvetica');
    lcnt = lcnt - br * 0.05;
end

end
