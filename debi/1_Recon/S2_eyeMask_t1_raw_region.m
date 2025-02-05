clear all; clc;
addpath(genpath('/media/debi/1,0 TB Disk/Backup/Recon_fork'));

%%

subject_num = 3;
datasetDir = '/home/debi/jaime/acquisitions/MREyeTrack/';
reconDir = '/media/debi/1,0 TB Disk/241120_JB/';
ETDir = '/home/debi/jaime/acquisitions/MREyeTrack/masks/';

% 0:up 1:down 2:left 3:right 4:center mask

% 'subject_1_mask_clean_0.3_0.mat'
if subject_num == 1
    datasetDir = [datasetDir, 'MREyeTrack_subj1/RawData_MREyeTrack_Subj1/'];  
elseif subject_num == 2
    datasetDir = [datasetDir, 'MREyeTrack_subj2/RawData_MREyeTrack_Subj2/'];
elseif subject_num == 3
    datasetDir = [datasetDir, 'MREyeTrack_subj3/RawData_MREyeTrack_Subj3/'];   
else
    datasetDir = [datasetDir, ' ']; 
end


if subject_num == 1
    otherDir = [reconDir, '/Sub001/T1_LIBRE_Binning/other/'];
elseif subject_num == 2
    otherDir = [reconDir, '/Sub002/T1_LIBRE_Binning/other/'];
elseif subject_num == 3
    otherDir = [reconDir, '/Sub003/T1_LIBRE_Binning/other/'];
else
    otherDir = [reconDir, '/Sub004/T1_LIBRE_Binning/other/'];
end

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

%%

for region_idx = 0:4
    switch region_idx
        case 0 
            disp('up mask')
        case 1 
            disp('down mask')
        case 2 
            disp('left mask')
        case 3 
            disp('right mask')
        case 4 
            disp('center mask')
    end


    th_ratio = 3/4;
    
    nShotOff = 15; 
    nSeg = 22; 
    eMask = eyeGenerateBinning(datasetDir,nShotOff, nSeg, th_ratio, ETDir, true);
    
    % Saving data and Convert to Monalisa format
    %--------------------------------------------------------------------------    
    
    %save(fullfile(param.savedir,'cMask.mat'),'param');
    %disp(['param is saved here:', param.savedir, '\cMask.mat'])
    
    eMaskFilePath = [otherDir,sprintf('eMask_th%.2f_raw_region%d.mat', th_ratio, region_idx)];
    
    % Save the CMask to the .mat file
    save(eMaskFilePath, 'eMask');
    disp('eMask has been saved here:')
    disp(eMaskFilePath)

end