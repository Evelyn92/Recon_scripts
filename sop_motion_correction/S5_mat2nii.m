%% Author: Yiwei
% This function converts mat file into nifti nii format
% for the following control by spm
% Feb 23. 2025
% Define the path to your .mat file and the output NIfTI file
%% Drag the .mat file into matlab 
% and save x into volume_data

volume_data = x;

%% Replace with your desired output NIfTI file path
nifti_folder = fullfile("/home/debi/yiwei/recon_results/250127_recon/", ...
    "Sub001/T1_LIBRE_Binning/output/mask_5s", ...
    'nii');       

if isfolder(nifti_folder)
    disp('The nifti_folder will be: ')
    disp(nifti_folder)
else
    mkdir(nifti_folder);
    disp('The nifti_folder is created: ')
    disp(nifti_folder)
end

nvolume = size(volume_data, 1);
disp(['number of volumes: ',num2str(nvolume)]);
%%
for idx = 1:nvolume
    vol_i = volume_data{idx};
    % Define NIfTI metadata (optional but recommended for completeness)
    % You can adjust these properties according to your needs.

    % NifTi file for each volume
    nifti_file = fullfile(nifti_folder, strcat('volume_', num2str(idx), '.nii'));
    niftiwrite(vol_i, nifti_file);
    disp(['NIfTI Data has been saved in folder: ', nifti_file]);
end



