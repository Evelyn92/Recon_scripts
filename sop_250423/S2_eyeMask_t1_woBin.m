% =====================================================
% Author: Yiwei Jia
% Date: April 17
% ------------------------------------------------
% Update:
% =====================================================
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));

% meas_MID00240_FID14889_sub5_main.dat    meas_MID00253_FID14902_sub2_pre_BC.dat
% meas_MID00241_FID14890_sub5_pre_HC.dat  meas_MID00254_FID14903_sub3_main.dat
% meas_MID00250_FID14899_sub5_pre_BC.dat  meas_MID00255_FID14904_sub4_main.dat
% meas_MID00251_FID14900_sub2_main.dat    meas_MID00256_FID14905_sub4_pre_HC.dat
% meas_MID00252_FID14901_sub2_pre_HC.dat  meas_MID00257_FID14906_sub4_pre_BC.dat
%%
subject_num = 1;

datasetDir = '/home/debi/yiwei/mreye_dataset/250423/';
reconDir = '/home/debi/yiwei/recon_results/250423/';
nShotOff = 14; 
nShot = 2055;
nSeg = 22; 

otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_woBinning/other/'];

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end



eMask = ones(1, nShot*nSeg);
eMask = (eMask>0);
eMask(1:nShotOff*nSeg) = 0;
% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    
eMaskFilePath = [otherDir,'eMask_woBin.mat'];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)

%%

