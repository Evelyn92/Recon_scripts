function [img_norm] = norm_single_volume(img)
img = abs(img);
img_norm = (img-min(img(:)))/(max(img(:))-min(img(:)));
disp('norm single volume done!')
end

