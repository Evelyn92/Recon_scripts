% =====================================================
% Author: Yiwei Jia
% Date: June 23
% ------------------------------------------------
% This script is used for generate binning mask 
% according to the ET mask, where
% sampling rate of ET mask: 1ms
% sampling rate of readouts: TR=8ms
% =====================================================
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));


%%
subject_num = 1;
datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-01/';
reconDir = '/home/debi/jaime/tmp/250613_JB/';
otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_Binning/other/'];

ETDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/masks_1206/';

% 0:up 1:down 2:left 3:right 4:center mask




% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

%%
% nShotOff should be aligned with the case of woBinning 
nShotOff = 14; 
nSeg = 22;

% winLen: the length of the readout sliding window to determine the preservation.
% th_ratio: the ratio for thresholding the ET mask.
winLen = 10;
th_ratio = 3/4;

% To be specific:
% To determine if the current readout should be preserved or not,
% we set the sliding window of winLen=10, so that we can check the
% corresponding ET window of length winLen*int(TR) = 80 ET points
% if more than th_ratio*winLen*int(TR) = 60 ET points are true (i.e. located within the criterion region)
    % we will maintain the current readout, i.e. binningMask value = 1
% else
    % discard this readout, i.e. binningMask value = 0
%%
% This function will guide you manually select the raw data and ET masks
% to generate the ET-guided binning mask for monalisa recon
% If you'd like to generate 4 bins with 4 masks,
% please enter nBin=4, and select ET masks for 4 times
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

