clear all; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_fork/'));

%%
subject_num=1;

datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot';
reconDir = '/home/debi/jaime/tmp/250613_JB';

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

if subject_num == 1
    bodyCoilFile    = [datasetDir, '/meas_MID2400614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile   = [datasetDir, '/meas_MID00615_FID182869_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile     = [datasetDir, '/meas_MID00605_FID182859_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 2
    bodyCoilFile    = [datasetDir, '/meas_MID00589_FID182843_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile   = [datasetDir, '/meas_MID00590_FID182844_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile     = [datasetDir, '/meas_MID00580_FID182834_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 3
    bodyCoilFile    = [datasetDir, '/meas_MID00563_FID182817_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile   = [datasetDir, '/meas_MID00564_FID182818_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile     = [datasetDir, '/meas_MID00554_FID182808_BEAT_LIBREon_eye_(23_09_24).dat'];
end

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
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
nShotOff = 14; 
nSeg = 22; 
% eMask = eyeGenerateBinning(datasetDir,nShotOff, nSeg, th_ratio, ETDir, true);
winLen = 10;
eMask = eyeGenerateBinningWin(datasetDir, nShotOff, nSeg,th_ratio, ETDir, winLen, true);
% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    

%save(fullfile(param.savedir,'cMask.mat'),'param');
%disp(['param is saved here:', param.savedir, '\cMask.mat'])

eMaskFilePath = [otherDir, sprintf('eMask_th%.2f_raw.mat', th_ratio)];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)

