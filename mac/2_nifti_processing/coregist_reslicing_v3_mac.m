addpath('/Users/cag/Documents/forclone/spm')
%%
subject_num = 4;
% Define image paths
switch subject_num
    case 1
        nii_folder = '/Users/cag/Documents/Dataset/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/NIFTI/';
        % Libre image
        a_image_name = 'DICOM_BEAT_LIBREon_eye_20230928191441_8.nii';
        % Vibe image
        b_image_name = 'DICOM_t1_vibe_fs_tra_p2_0.4_HR_64_canaux_20230928191441_9.nii'; 
        % Mprage image
        c_image_name = 'DICOM_wip19_mprage_1iso_cs4p2_20230928191441_5.nii'; 
    case 2    
        nii_folder = '/Users/cag/Documents/Dataset/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/MR_EYE_Subj02/NIFTI/';
        % Libre image
        a_image_name = 'DICOM_BEAT_LIBREon_eye_20230926164342_8001.nii'; 
        % Vibe image
        b_image_name = 'DICOM_t1_vibe_fs_tra_p2_0.4_HR_64_canaux_20230926164342_9001.nii'; 
        % Mprage image
        c_image_name = 'DICOM_wip19_mprage_1iso_cs4p2_20230926164342_5001.nii'; 
    case 4
        nii_folder = '/Users/cag/Documents/Dataset/1_Pilot_MREye_Data/Sub004/230923_anatomical_MREYE_study/MR_EYE_Subj04/NIFTI/';
        % Libre image
        a_image_name = 'DICOM_BEAT_LIBREon_eye_20230928151925_8.nii'; 
        % Vibe image
        b_image_name = 'DICOM_t1_vibe_fs_tra_p2_0.4_HR_64_canaux_20230928151925_9.nii'; 
        % Mprage image
        c_image_name = 'DICOM_wip19_mprage_1iso_cs4p2_20230928151925_5.nii'; 
   
    
end

a_image = fullfile(nii_folder, a_image_name);
b_image = fullfile(nii_folder, b_image_name);
c_image = fullfile(nii_folder, c_image_name);
%% Inspect the header before co-registration
V_a_1 = spm_vol(a_image);
disp('V_a_1'); disp(V_a_1);
V_b_1 = spm_vol(b_image);
disp('V_b_1'); disp(V_b_1);

% Co_register
co_register(b_image, a_image)

% Inspect the header after co-registration
V_a_2 = spm_vol(a_image);
disp('V_a_2'); disp(V_a_2);
V_b_2 = spm_vol(b_image);
disp('V_b_2'); disp(V_b_2)

disp(['V_b remains the same? ', num2str(isequal(V_b_1.mat, V_b_2.mat))])
%%
%--------------------------------------------------------------------------
% Reslicing part
%--------------------------------------------------------------------------
reslicing(b_image, a_image)

