clc;
clear all;
%% Inspect LIBRE W/WO Binning
x_Binning = load('/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub002/T1_LIBRE_Binning/output/th7/x_nIter3.mat');
% /media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub001/T1_LIBRE_woBinning/output/x_nIter3.mat
% x_Binning = x;
% clear x

x_woBinning = load('/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub002/T1_LIBRE_woBinning/output/x_nIter3.mat');


if ~ isequal(x_Binning, x_woBinning)
    disp('The recon results w/wo binning is not equal!')
end

disp('Slice the two results.')
x_Binning = x_Binning.x;
x_woBinning = x_woBinning.x;
slice_wobinning = x_woBinning{1};
slice_binning = x_Binning{1};

close all;
% set the index of slice for inspecting
% xLeftEyeRange = 60:120;
% yLeftEyeRange = 130:180;
xLeftEyeRange = (240-120):(240-60);
yLeftEyeRange = (240-180):(240-130);
z=112;
bmImage(slice_binning(:,:,:))
bmImage(slice_binning(xLeftEyeRange,yLeftEyeRange,:))
%%
bmImage(slice_wobinning(:,:,:))
bmImage(slice_wobinning(xLeftEyeRange,yLeftEyeRange,:))
%%
% close all;
% bmImage(slice_binning(xLeftEyeRange,yLeftEyeRange,z)-slice_wobinning(xLeftEyeRange,yLeftEyeRange,z))

%%
slice_flipped_binning = flip(slice_binning);
slice_flipped_wobinning = flip(slice_wobinning);

%%
close all;
% set the index of slice for inspecting
% xLeftEyeRange = 60:120;
% yLeftEyeRange = 130:180;
xLeftEyeRange = (60):(120);
yLeftEyeRange = (130):(180);
z=112;
bmImage(slice_flipped_binning(:,:,:))
bmImage(slice_flipped_wobinning(xLeftEyeRange,yLeftEyeRange,:))
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
