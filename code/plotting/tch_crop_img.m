function cropped_img = tch_crop_img(img, bg_color)

if nargin < 2; bg_color = [255 255 255]; end
if length(bg_color) > 3; error('Invalid background color'); end
if size(img, 3) > 3; error('Image format not recognized'); end
if size(img, 3) < 3; img = repmat(img(:, :, 1), 1, 1, 3); end
[img_h, img_w, ~] = size(img); bg_color = reshape(bg_color, 1, 1, 3);
bcol = repmat(bg_color, img_h, 1, 1); brow = repmat(bg_color, 1, img_w, 1);

rcnt_t = 1; rcnt_b = img_h; ccnt_l = 1; ccnt_r = img_w;
while isequal(img(rcnt_t, :, :), brow); rcnt_t = rcnt_t + 1; end
while isequal(img(rcnt_b, :, :), brow); rcnt_b = rcnt_b - 1; end
while isequal(img(:, ccnt_l, :), bcol); ccnt_l = ccnt_l + 1; end
while isequal(img(:, ccnt_r, :), bcol); ccnt_r = ccnt_r - 1; end
cropped_img = img(rcnt_t:rcnt_b, ccnt_l:ccnt_r, :);

end
