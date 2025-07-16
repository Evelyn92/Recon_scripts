%%
subject_num = 2;
% Define image paths

if subject_num == 1
    nii_folder = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/NIFTI_NII_origin/';
    png_folder = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/png/t1_vibe_origin/';
    % Libre image
    a_image_name = 'DICOM_BEAT_LIBREon_eye_20230928191441_8.nii';
    % Vibe image
    b_image_name = 'resliced_DICOM_t1_vibe_fs_tra_p2_0.4_HR_64_canaux_20230928191441_9.nii'; 
else
    nii_folder = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/MR_EYE_Subj02/NIFTI_NII_origin/';
    png_folder = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/MR_EYE_Subj02/png/t1_vibe_origin/';
    % Libre image
    a_image_name = 'DICOM_BEAT_LIBREon_eye_20230926164342_8001.nii'; 
    % Vibe image
    b_image_name = 'resliced_DICOM_t1_vibe_fs_tra_p2_0.4_HR_64_canaux_20230926164342_9001.nii'; 
    
end

a_image = fullfile(nii_folder, a_image_name);
b_image = fullfile(nii_folder, b_image_name);

%%
% Load NIfTI file
V = spm_vol(b_image);

% Read voxel data
imageData = spm_read_vols(V);

% Display basic information
disp(['Image Dimensions: ', num2str(size(imageData))]);

%%
% Specify the output directory

if ~exist(png_folder, 'dir')
    mkdir(png_folder);
end

% Number of slices
numSlices = size(imageData, 3);

% Save each slice as PNG
for slice = 1:numSlices
    % Extract the slice
    sliceData = imageData(:,:,slice);
    
    % Normalize for display
    sliceData = mat2gray(sliceData); % Convert to grayscale [0, 1]
    rotatedData = rot90(sliceData, -1); % -1 means 90 degrees clockwise
    % Save as PNG
    pngFilename = fullfile(png_folder, sprintf('slice_%03d.png', slice));
    imwrite(rotatedData, pngFilename);
end

disp('All slices have been saved as PNG files.');



