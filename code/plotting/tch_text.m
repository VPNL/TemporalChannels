function txt = tch_text(x, y, str, fs, col, h_align, v_align, fw)

if nargin < 3; error('Not enough input arguments.'); end
if nargin < 4; fs = 7; end
if nargin < 5; col = [0 0 0]; end
if nargin < 6; h_align = 'c'; end
if nargin < 7; v_align = 'm'; end
if nargin < 8; fw = 'normal'; end

switch h_align(1)
    case 'l'
        ha = 'left';
    case 'r'
        ha = 'right';
    otherwise
        ha = 'center';
end
switch v_align(1)
    case 't'
        va = 'top';
    case 'b'
        va = 'bottom';
    otherwise
        va = 'middle';
end

txt = text(x, y, str, 'FontName', 'Helvetica', 'FontSize', fs, 'Color', col, ...
    'HorizontalAlignment', ha, 'VerticalAlignment', va, 'FontWeight', fw);


end
