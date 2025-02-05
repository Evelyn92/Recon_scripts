clear all; clc;
addpath(genpath('/media/debi/1,0 TB Disk/Backup/Recon_fork'));

%%
datasetDir = '/home/debi/jaime/acquisitions/MREyeTrack/';
reconDir = '/media/debi/1,0 TB Disk/241120_JB/';

for subject_num = 2
if subject_num == 1
    datasetDir = [datasetDir, 'MREyeTrack_subj1/RawData_MREyeTrack_Subj1/'];
    ETDir = [datasetDir, ' '];
elseif subject_num == 2
    datasetDir = [datasetDir, 'MREyeTrack_subj2/RawData_MREyeTrack_Subj2/'];
    ETDir = [datasetDir, ' '];
elseif subject_num == 3
    datasetDir = [datasetDir, 'MREyeTrack_subj3/RawData_MREyeTrack_Subj3/'];
    ETDir = [datasetDir, ' '];
else
    datasetDir = [datasetDir, ' '];
    ETDir = [datasetDir, ' '];
end


if subject_num == 1
    bodyCoilFile     = [datasetDir, '/meas_MID2400614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00615_FID182869_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00605_FID182859_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 2
    bodyCoilFile     = [datasetDir, '/meas_MID00589_FID182843_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00590_FID182844_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00580_FID182834_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 3
    bodyCoilFile = [datasetDir, '/meas_MID00563_FID182817_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile = [datasetDir, '/meas_MID00564_FID182818_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00554_FID182808_BEAT_LIBREon_eye_(23_09_24).dat'];
else
    % bodyCoilFile = [datasetDir, '/meas_MID00349_FID57815_BEAT_LIBREon_eye_BC_BC.dat'];
    % arrayCoilFile = [datasetDir, '/meas_MID00350_FID57816_BEAT_LIBREon_eye_HC_BC.dat'];
    % measureFile = [datasetDir, '/meas_MID00333_FID57799_BEAT_LIBREon_eye.dat'];
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


nShotOff = 15; 
nSeg = 22; 
eMask = eyeGenerateWoBinning(datasetDir, nShotOff, nSeg);

% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    

%save(fullfile(param.savedir,'cMask.mat'),'param');
%disp(['param is saved here:', param.savedir, '\cMask.mat'])

eMaskFilePath = [otherDir,'eMask_woBin.mat'];

% Save the CMask to the .mat file
save(eMaskFilePath, 'eMask');
disp('eMask has been saved here:')
disp(eMaskFilePath)


end
%%

