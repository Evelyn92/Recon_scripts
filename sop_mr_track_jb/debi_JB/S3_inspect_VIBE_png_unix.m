clc;
clear all;
%% Inspect VIBE image

%%
sub001_t1vibe_png = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/png/t1_vibe';
k = 114;
img_k = imread(fullfile(sub001_t1vibe_png, sprintf('slice_%03d.png', k)));
img_k_rot = imrotate(img_k, -90);
img_k_scale = imresize(img_k_rot, 0.8);
img_k_scale = mat2gray(img_k_scale); % Normalizes the image
% imshow(img_k_rot);
% imshow(img_k_scale);
imshow(img_k_scale);
% up:down, left:right
