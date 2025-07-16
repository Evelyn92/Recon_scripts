% >> load('/home/debi/yiwei/recon_results/250603/Sub001/T1_LIBRE_woBinning/output/mask_idea_jb_rect_pulse/xrms.mat')
% >> xrms_i = permute(xrms, [2,3,1]);
% >> load('/home/debi/yiwei/recon_results/250603/Sub002/T1_LIBRE_woBinning/output/mask_pulseq_rect_pulse/xrms.mat')
% >> xrms_p = permute(xrms, [2,3,1]);
% >> bmImage(xrms_p)
% >> xrms_i_p = cat(2,[xrms_i, xrms_p]);
% >> bmImage(xrms_i_p);
% Ensure the images are in double precision

x1 = permute(x_libre_6p2, [1,3,2]);
x1 = flip(x1,1);
bmImage(x1)

x2 = permute(x_gre_6p2, [1,3,2]);
x2 = flip(x2,1);
bmImage(x2)
%%
[img1, img2] = norm_two_image(x1(61:120,:,:),x2(61:120,:,:));
img_1_2 = cat(1,img1, img2);
bmImage(img_1_2);

%% sl_start = 40;
show_image = img_1_2(:,:,sl_start);
for slice = 45:5:80
    show_image = cat(2,[show_image, img_1_2_vertical(:,:,slice)]);
end
bmImage(show_image);

%%
diff_map = diff_volume(xrms_2,xrms_5);
sl_start = 100;
show_diff_image = permute(diff_map(:,:,sl_start), [2,1,3]);
for slice_idx = 110:10:140
show_diff_image = cat(2,[show_diff_image, permute(diff_map(:,:,slice_idx), [2,1,3])]);

end
loLev = 0.0;upLev = 0.6;[imClip, rgb_vec] = relaxationColorMap('T1', show_diff_image, loLev, upLev);
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