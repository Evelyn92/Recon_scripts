clear all; clc;
% addpath(genpath("/media/sinf/1,0 TB Disk/Backup/Recon_fork"));
addpath(genpath('/home/debi/yiwei/forclone/Recon_fork/'));
%%
% Coil sensitivity -> binning mask cMask -> Mitosius
% ------------------------------------------------------------------
% ------------------------------------------------------------------

%% Initialize the directories and acquire the Coil
subject_num = 1;
%
datasetDir = '/home/debi/jaime/mreye_track/data/pilot/sub-01/';
reconDir = '/home/debi/jaime/tmp/250613_JB/';

if subject_num == 1
    datasetDir = [datasetDir, 'rawdata'];
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
    bodyCoilFile     = [datasetDir, '/meas_MID00614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
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
%% Start coil estimation
disp('-------------------------------------')
disp('Start Coil Estimation')
disp('-------------------------------------')
% Read metadata from the twix
bmTwix_info(bodyCoilFile)
bmTwix_info(arrayCoilFile)
%%
inspect_raw_data = true;
if inspect_raw_data
    bmTwix_info(measureFile);
    check_data_timestamp(datasetDir,false);
end

%% The input initialization is not automatizable, since the name/content of
% the variables depends on the sequence programmed. For the moment we just
% read and print the twixinfo and you need to adjut the parameters yourself

% Maybe: write a function for automate ISMR rawData format(Standard) reading.
% All trajectory information, to generate the trajectory. 
N            = 480; 
nSeg         = 22; 
nShot        = 419; 
% The FoV can be set as 48 here (low-resolution), 
% since we want to estimate coil sensitivity

acquisitonFoV          = [240, 240, 240];
% This number (nShotOff) has to be adapted based on the observation of the 
% steady-state graph
nShotOff     = 15; 
N_u          = [48, 48, 48]; % the size of k-space
reconFov = 240;
dK_u         = [1, 1, 1]./reconFov; % interval between k-space

if subject_num == 1
    nCh_array    = 44;
elseif subject_num == 2
    nCh_array    = 52;
elseif subject_num == 3
    nCh_array = 52;
else
    nCh_array = 999 ;
end

nCh_body     = 2; 
% We will need to ask for a predefined format. Hence we need to read the data outside and
% pass the y_body, t has to be computed by the user.  
% We use trajectory in using the phisical dimentions without convention [-0.5,0.5]
% You have to see how much data needs to be discarded by looking at the
% graph, you don't want to keep data if you are not in steady-state.
%%
[y_body, t, ve] = bmCoilSense_nonCart_dataFromTwix( bodyCoilFile, ...
                                                    N_u, ...
                                                    N, ...
                                                    nSeg, ...
                                                    nShot, ...
                                                    nCh_body, ...
                                                    reconFov, ...
                                                    nShotOff);
% same for arraycoil 
y_array         = bmCoilSense_nonCart_dataFromTwix( arrayCoilFile, ...
                                                    N_u, ...
                                                    N, ...
                                                    nSeg, ...
                                                    nShot, ...
                                                    nCh_array, ...
                                                    reconFov, ...
                                                    nShotOff);
% compute the gridding matrices (Gn = approximation of inverse, Gu = Forward,
% Gut = transposed of Gu) Gn and Gut are both backward
[Gn, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u); 
disp('Gn, Gu, Gut are computed')
%
% You need to Reassign the xmin, xmax & ymin, ymax & zmin, zmax 
% To do it you need to run the function below (bmCoilSense_nonCart_mask)
% control + E: to change the tresholds
% shift + E: to set the constast chosen

%% Box excluding coordinates
x_min = 1; 
x_max = 40;

y_min = 9; 
y_max = 41;

z_min = 4; 
z_max = 48;

% Two thresholds
th_RMS = 19; 
th_MIP = 16; 

close_size = []; 
open_size  = []; 
% Mask computation, the two thresholds to exclude artifacts from the region
% where there is no signal, like air in the lungs: 
% 1 for the Root mean square
% 1 for the Maximum intensity projection

% Make some work to try to automate it. Define two scripts.
% (Automatic and Advanced)
m = bmCoilSense_nonCart_mask(   y_body, Gn, ...
                                x_min, x_max, ...
                                y_min, y_max, ...
                                z_min, z_max, ...
                                th_RMS, th_MIP, ...
                                close_size, ...
                                open_size, ...
                                true);
% Reset the parameters in the cell above and rerun this cell.
% Figure 1:
    % X axis: the normalized intensity 
    % Y axis: the percentage of pixels whose internsities are above the normalized intensity
    % The thresholds should be at the elbow of the curve (most of the pixels)

%%
close all;
% Select one body coil and compute its sensitivity
[y_ref, C_ref] = bmCoilSense_nonCart_ref(y_body, Gn, m, []); 

% Estimate the coil sensitivity of each surface coil using one body coil
% image as reference image C_c = (X_c./x_ref)
C_array_prime = bmCoilSense_nonCart_primary(y_array, y_ref, C_ref, Gn, ve, m);

% Do a recon, predending the selected body coil is one channel among the
% others, and optimize the coil sensitivity estimate by alternating steps
% Of gradient descent (X,C)
nIter = 5; 
[C, x] = bmCoilSense_nonCart_secondary(y_array, C_array_prime, y_ref, C_ref, Gn, Gu, Gut, ve, nIter, false); 
close all;
%% Save C into the folder
if subject_num == 1
    saveCDirList = {'/Sub001/T1_LIBRE_Binning/C/','/Sub001/T1_LIBRE_woBinning/C/'};
elseif subject_num == 2
    saveCDirList = {'/Sub002/T1_LIBRE_Binning/C/','/Sub002/T1_LIBRE_woBinning/C/'};
elseif subject_num == 3
    saveCDirList = {'/Sub003/T1_LIBRE_Binning/C/','/Sub003/T1_LIBRE_woBinning/C/'};
else
    saveCDirList = {'/Sub004/T1_LIBRE_Binning/C/','/Sub004/T1_LIBRE_woBinning/C/'};
end

for idx = 1:2
    saveCDir     = [reconDir, saveCDirList{idx}];
    CfileName = 'C.mat';
    
    % Create the folder if it doesn't exist
    if ~exist(saveCDir, 'dir')
        mkdir(saveCDir);
    end
    
    % Full path to  C file
    CfilePath = fullfile(saveCDir, CfileName);
    
    % Save the matrix C to the .mat file
    save(CfilePath, 'C');
    disp('Coil sensitivity C has been saved here:')
    disp(CfilePath)
end




