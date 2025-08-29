% >> load('/home/debi/yiwei/recon_results/250603/Sub001/T1_LIBRE_woBinning/output/mask_idea_jb_rect_pulse/xrms.mat')
% >> xrms_i = permute(xrms, [2,3,1]);
% >> load('/home/debi/yiwei/recon_results/250603/Sub002/T1_LIBRE_woBinning/output/mask_pulseq_rect_pulse/xrms.mat')
% >> xrms_p = permute(xrms, [2,3,1]);
% >> bmImage(xrms_p)
% >> xrms_i_p = cat(2,[xrms_i, xrms_p]);
% >> bmImage(xrms_i_p);
% Ensure the images are in double precision
bmImage(xlibre3);
bmImage(xlibre4);
%%
x1 = permute(xlibre3, [1,3,2]);
x1 = flip(x1,1);
x1_rot = rotate_volume(x1, -12);
bmImage(x1_rot)

x2 = permute(xlibre4, [1,3,2]);
x2 = flip(x2,1);
x2_rot = rotate_volume(x2, -12);
bmImage(x2_rot)
%%
x3 = permute(xgre3, [1,3,2]);
x3 = flip(x3,1);
x3_rot = rotate_volume(x3, -12);
bmImage(x3_rot)
% 
x4 = permute(xgre4, [1,3,2]);
x4 = flip(x4,1);
x4_rot = rotate_volume(x4, -12);
bmImage(x4_rot)
%%
[img1, img2] = norm_two_image(x3_rot(range{1}*2, range{2}*2, range{3}*2),x4_rot(range{1}*2, range{2}*2, range{3}*2));
img_1_2 = cat(1,img1, img2);
bmImage(img_1_2);

%% 
sl_start = 38;
inc=4;
sl_end=54;

show_image = img_1_2(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2(:,:,slice)]);
end
bmImage(show_image);

%%
diff_map = diff_volume(img1,img2);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);

end
loLev = 0.05;upLev = 0.4;[imClip, rgb_vec] = relaxationColorMap('T1', show_diff_image, loLev, upLev);
figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(imClip, 'DisplayRange', [loLev, upLev], 'InitialMagnification', 'fit'); 
% title(strcat('difference map-', num2str(slice_idx)));
colormap(rgb_vec); colorbar;
%%
function [img1_scaled, img2_scaled] = norm_two_image(img1,img2)
img1 = double(abs(img1));
img2 = double(abs(img2));

% Scale img1 to [0, 1]
img1_scaled = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));

% Scale img2 to [0, 1]
img2_scaled = (img2 - min(img2(:))) / (max(img2(:)) - min(img2(:)));
disp('Scaling Done')
end

function [diff_img]=diff_volume(img1,img2)
[img1_s,img2_s] = norm_two_image(img1,img2);
diff_img = abs(img1_s-img2_s);
disp(max(diff_img(:)))
end

function V_rot = rotate_volume(V,theta)
V = double(abs(V));

% theta: degree
% positive: counterclockwise; negative: clockwise
    V_rot = zeros(size(V));
    for k = 1:size(V,3)
        V_rot(:,:,k) = imrotate(V(:,:,k), theta, 'bilinear', 'crop');  % or 'loose'
    end
end