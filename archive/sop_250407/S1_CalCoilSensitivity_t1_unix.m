% =====================================================
% Author: Yiwei Jia
% Date: Feb 5
% ------------------------------------------------
% [Coil sensitivity] -> binning mask eMask -> Mitosius
% Update: this script is derived from Demo script
% by Mauro in Monalisa version Feb.5
% The old script has issue when running mask generation
% With readers, the param setting is more organized
% =====================================================
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));



%% Initialize the directories and acquire the Coil
subject_num = 1;
%
datasetDir = '/home/debi/yiwei/mreye_dataset/250127_acquisition/';
reconDir = '/home/debi/yiwei/recon_results/250127_recon/';
if subject_num == 1
    datasetDir = [datasetDir, ''];
    ETDir = [datasetDir,'ET_EDF/'];
elseif subject_num == 2
    datasetDir = [datasetDir, ' '];
    ETDir = [datasetDir, ' '];
elseif subject_num == 3
    datasetDir = [datasetDir, ' '];
    ETDir = [datasetDir, ' '];
else
    datasetDir = [datasetDir, ' '];
    ETDir = [datasetDir, ' '];
end

if subject_num == 1
    bodyCoilFile     = [datasetDir, '/meas_MID00345_FID214641_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00346_FID214642_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00332_FID214628_BEAT_LIBREon_eye_(23_09_24)_sc_trigger.dat'];
elseif subject_num == 2
    bodyCoilFile     = [datasetDir, '..'];
    arrayCoilFile    = [datasetDir, '..'];
    measureFile = [datasetDir, '..'];
elseif subject_num == 3
    bodyCoilFile = [datasetDir, '/.dat'];
    arrayCoilFile = [datasetDir, '/.dat'];
    measureFile = [datasetDir, '/.dat'];
else
    % bodyCoilFile = [datasetDir, '/meas_MID00349_FID57815_BEAT_LIBREon_eye_BC_BC.dat'];
    % arrayCoilFile = [datasetDir, '/meas_MID00350_FID57816_BEAT_LIBREon_eye_HC_BC.dat'];
    % measureFile = [datasetDir, '/meas_MID00333_FID57799_BEAT_LIBREon_eye.dat'];
end


%% Load and Configure Data
% Read data using the library's `createRawDataReader` function
% This readers makes the usage of Siemens and ISMRMRD files equivalent for
% the libraray
bodyCoilreader = createRawDataReader(bodyCoilFile, false);
arrayCoilReader = createRawDataReader(arrayCoilFile, false);

% Ensure consistency in number o1f shot-off points
nShotOff = max(bodyCoilreader.acquisitionParams.nShot_off, ...
               arrayCoilReader.acquisitionParams.nShot_off);
bodyCoilreader.acquisitionParams.nShot_off = nShotOff;
arrayCoilReader.acquisitionParams.nShot_off = nShotOff;

%% Parameters
dK_u = [1, 1, 1] ./ arrayCoilReader.acquisitionParams.FoV;   % Cartesian grid spacing
N_u = [48, 48, 48];             % Adjust this value as needed
%% Compute Trajectory and Volume Elements
[y_body, t, ve] = bmCoilSense_nonCart_data(bodyCoilreader, N_u);
y_surface = bmCoilSense_nonCart_data(arrayCoilReader, N_u);

% Compute the gridding matrices (subscript is a reminder of the result)
% Gn is from uniform to Non-uniform
% Gu is from non-uniform to Uniform
% Gut is Gu transposed
[Gn, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);
%% Create Mask
mask = bmCoilSense_nonCart_mask_automatic(y_body, Gn, false);
%% Estimate Coil Sensitivity
% Reference coil sensitivity using the body coils. This is used as 
% a reference to estiamte the sensitivity of each head coil
[y_ref, C_ref] = bmCoilSense_nonCart_ref(y_body, Gn, mask, []);

% Head coil sensitivity estimate using body coil reference
C_array_prime = bmCoilSense_nonCart_primary(y_surface, y_ref, C_ref, Gn, ve, mask);
%% Refine the sensitivity estimate with optimization
nIter = 5;
[C, x] = bmCoilSense_nonCart_secondary(y_surface, C_array_prime, y_ref, ...
                                       C_ref, Gn, Gu, Gut, ve, nIter, false);

%% Display Results
bmImage(C);
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




