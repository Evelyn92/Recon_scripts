

%% only apply for 0829
x_fa = cat(2, [x1,x2,x4]); %compare FA 4 8 16
x_freq_sign = cat(2, [x1,x3]); %compare freq offset neg pos
x_swap_10 = cat(2, [x1,x5]); %compare grad swap or not

%%
close all;clc
x = {x1, x2, x3, x4, flip(flip(flip(x5,1),2),3), x6};

%%

x1_trans = norm_image(permute(x{1}, [2,1,3]));
x1_sag = norm_image(flip(rot90(permute(x1_trans, [1,3,2]), 1),1));
x1_coronal = norm_image(permute(x1_trans, [3,2,1]));
% bmImage(x1_trans);
% bmImage(x1_sag);
% bmImage(x1_coronal);
%
x2_trans = norm_image(permute(x{6}, [2,1,3]));
x2_sag = norm_image(flip(rot90(permute(x2_trans, [1,3,2]),1),1));
x2_coronal = norm_image(permute(x2_trans, [3,2,1]));
% bmImage(x2_trans);
% bmImage(x2_sag);
% bmImage(x2_coronal);
%
% x3_trans = norm_image(permute(x{4}, [2,1,3]));
% x3_sag = norm_image(flip(rot90(permute(x3_trans, [1,3,2]),1),1));
% x3_coronal = norm_image(permute(x3_trans, [3,2,1]));
% bmImage(x3_trans);
% bmImage(x3_sag);
% bmImage(x3_coronal);

%
[img1_trans, img2_trans] = norm_two_image(x1_trans,x2_trans);
img_1_2_trans = cat(1,img1_trans, img2_trans);
bmImage(img_1_2_trans)

[img1_sag, img2_sag] = norm_two_image(x1_sag,x2_sag);
img_1_2_sag = cat(1,img1_sag, img2_sag);
bmImage(img_1_2_sag)

% img_trans_sag = cat(2,img_1_2_trans, img_1_2_sag);
% bmImage(img_trans_sag);

%% ===== transverse ====================
close all;
sl_start = 50;
inc=3;
sl_end=70;
% show_image = img1_trans(:,:,sl_start);
show_image = img_1_2_trans(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2_trans(:,:,slice)]);
end
bmImage(show_image);

%%
diff_map = diff_volume(img1_trans,img2_trans);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);
end

figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(show_diff_image);
colorbar;colormap('redblue'); caxis([-max(abs(diff_map(:))), max(abs(diff_map(:)))]);

% loLev = -1;upLev = 0.1;[imClip, rgb_vec] = relaxationColorMap('T1', show_diff_image, loLev, upLev);
% figure('Color', 'white'); set(gca, 'Color', 'white'); 
% imshow(imClip, 'DisplayRange', [loLev, upLev], 'InitialMagnification', 'fit'); 
% % title(strcat('difference map-', num2str(slice_idx)));
% colormap(rgb_vec); colorbar;

%% ===== sagittal left====================

sl_start = 50;
inc=5;
sl_end=80;

show_image = img_1_2_sag(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2_sag(:,:,slice)]);
end
bmImage(show_image);

%%
diff_map = diff_volume(img1_sag,img2_sag);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);
end

figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(show_diff_image);
colorbar;colormap('redblue'); caxis([-max(abs(diff_map(:)))/4*3, max(abs(diff_map(:)))/4*3]);

%% ====coronal ====
sl_start = 94;
inc=4;
sl_end=110;

show_image = img_1_2_sag(:,:,sl_start);

for slice = (sl_start+inc):inc:sl_end
    show_image = cat(2,[show_image, img_1_2_sag(:,:,slice)]);
end
bmImage(show_image);


diff_map = diff_volume(img1_sag,img2_sag);
show_diff_image = diff_map(:,:,sl_start);
for slice_idx = sl_start+inc:inc:sl_end
show_diff_image = cat(2,[show_diff_image, diff_map(:,:,slice_idx)]);
end

figure('Color', 'white'); set(gca, 'Color', 'white'); 
imshow(show_diff_image);
colorbar;colormap('redblue'); caxis([-max(abs(diff_map(:)))/4*3, max(abs(diff_map(:)))/4*3]);


%%
function [img1_scaled] = norm_image(img1)
img1 = double(abs(img1));
% Scale img1 to [0, 1]
img1_scaled = (img1 - min(img1(:))) / (max(img1(:)) - min(img1(:)));

disp('Scaling Done')
end

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
diff_img = img1_s-img2_s;
disp(max(diff_img(:)))
disp(min(diff_img(:)))
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