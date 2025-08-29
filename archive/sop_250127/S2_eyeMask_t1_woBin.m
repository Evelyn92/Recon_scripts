% =====================================================
% Author: Yiwei Jia
% Date: Feb 5
% ------------------------------------------------
% Update:
% =====================================================
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));


%%
subject_num = 1;

datasetDir = '/home/debi/yiwei/mreye_dataset/250127_acquisition/';
reconDir = '/home/debi/yiwei/recon_results/250127_recon/';

if subject_num == 1
    datasetDir = [datasetDir, ''];
    ETDir = [datasetDir, ' '];
elseif subject_num == 2
    datasetDir = [datasetDir, ''];
    ETDir = [datasetDir, ' '];
elseif subject_num == 3
    datasetDir = [datasetDir, ''];
    ETDir = [datasetDir, ' '];
else
    datasetDir = [datasetDir, ' '];
    ETDir = [datasetDir, ' '];
end



if subject_num == 1
    otherDir = [reconDir, '/Sub001/T1_LIBRE_woBinning/other/'];
elseif subject_num == 2
    otherDir = [reconDir, '/Sub002/T1_LIBRE_woBinning/other/'];
elseif subject_num == 3
    otherDir = [reconDir, '/Sub003/T1_LIBRE_woBinning/other/'];
else
    otherDir = [reconDir, '/Sub004/T1_LIBRE_woBinning/other/'];
end

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end


nShotOff = 14; 
nSeg = 22; 
eMask = eyeGenerateWoBinning(datasetDir, nShotOff, nSeg);

% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    
eMaskFilePath = [otherDir,'eMask_woBin.mat'];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)



%%

