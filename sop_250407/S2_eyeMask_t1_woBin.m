% =====================================================
% Author: Yiwei Jia
% Date: April 7
% ------------------------------------------------
% Update:
% =====================================================
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));


%%
subject_num = 5;

datasetDir = '/home/debi/yiwei/mreye_dataset/250407/';
reconDir = '/home/debi/yiwei/recon_results/250407/';
nShotOff = 14; 
nShot = 419;
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

