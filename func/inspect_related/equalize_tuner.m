%left side
xLeftEyeRange = 120:240;
yLeftEyeRange = 1:120;
roi_range{1} = xLeftEyeRange;
roi_range{2} = yLeftEyeRange;
xRange = roi_range{1};
yRange = roi_range{2};
%%
nii_folder_1 = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/NIFTI_NII_origin/';
image_name_1 = 'resliced_DICOM_t1_vibe_fs_tra_p2_0.4_HR_64_canaux_20230928191441_9.nii';
sub001_t1vibe = [nii_folder_1,image_name_1];
%%
input_file = sub001_t1vibe;

recon_data = niftiread(input_file);
[~, ~, slices] = size(recon_data);
i=100;
slice = abs(recon_data(:, :, i));


rotated_slice = rot90(slice, -1);% -1 means 90 degrees clockwise
slice = abs(rotated_slice(xRange, yRange));


slice_norm = mat2gray(slice); % Normalizes the image
imshow(slice_norm)
figure, imhist(slice_norm);           % Show the histogram
%%
min_val = 0.01;  % Lower threshold
max_val = 0.8;  % Upper threshold

% Clip the intensity values
slice_I_clipped = slice_norm;
slice_I_clipped(slice_norm < min_val) = min_val;
slice_I_clipped(slice_norm > max_val) = max_val;

% Normalize to [0, 1] after clipping
slice_redistributed = (slice_I_clipped - min(slice_I_clipped(:))) / (max(slice_I_clipped(:)) - min(slice_I_clipped(:)));

% Display the results
figure;
subplot(1,2,1); imshow(slice_norm, []); title('Original Image');
subplot(1,2,2); imshow(slice_redistributed, []); title('Clipped & Redistributed Image');


figure;
subplot(1,2,1); imhist(slice_norm); title('Original Histogram');
subplot(1,2,2); imhist(slice_redistributed); title('Clipped & Redistributed Histogram');

