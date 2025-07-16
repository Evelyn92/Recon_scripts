clear all; clc;
addpath(genpath('/media/debi/1,0 TB Disk/Backup/Recon_fork'));

%%
subject_num=1;

datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot';
reconDir = '/home/debi/jaime/tmp/250613_JB';

% for the masks in ETDir --> 0:up 1:down 2:left 3:right 4:center mask
if subject_num == 1
    datasetDir = [datasetDir, '/sub-01/rawdata'];
    ETDir      = [datasetDir, '/masks_1206/Sub001'];
    otherDir   = [reconDir, '/Sub001/T1_LIBRE_Binning/other/'];
elseif subject_num == 2
    datasetDir = [datasetDir, '/sub-02/rawdata'];
    ETDir      = [datasetDir, '/masks_1206/Sub002'];
    otherDir   = [reconDir, '/Sub002/T1_LIBRE_Binning/other/'];
elseif subject_num == 3
    datasetDir = [datasetDir, '/sub-03/rawdata'];
    ETDir      = [datasetDir, '/masks_1206/Sub003'];
    otherDir   = [reconDir, '/Sub003/T1_LIBRE_Binning/other/'];
else
    datasetDir = [datasetDir, ' '];
    ETDir      = [datasetDir, ' '];
    otherDir   = [reconDir, ' '];
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
th_ratio = 3/4;
nShotOff = 15; 
nSeg = 22; 
% eMask = eyeGenerateBinning(datasetDir,nShotOff, nSeg, th_ratio, ETDir, true);
winLen = 10;
eMask = eyeGenerateBinningWin(datasetDir, nShotOff, nSeg,th_ratio, ETDir, winLen, true);
% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    

%save(fullfile(param.savedir,'cMask.mat'),'param');
%disp(['param is saved here:', param.savedir, '\cMask.mat'])

eMaskFilePath = [otherDir,sprintf('eMask_th%.2f_clean.mat', th_ratio)];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)

