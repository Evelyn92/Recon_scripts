% =====================================================
% Author: Yiwei Jia
% Date: March 20
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
subject_num = 4;
%
datasetDir = '/home/debi/yiwei/mreye_dataset/debug_0320/data/';
reconDir = '/home/debi/yiwei/recon_results/250320/';


if subject_num == 1
    %a: rf1: 2 rf pulses
    meas_name_suffix = '_MID0001_JB_LIBRE2p2_a8';
    hc_name_suffix = '_MID0001_JB_LIBRE2p2_a8';
    bc_name_suffix = '_MID0001_JB_LIBRE2p2_a8_BC';
elseif subject_num == 2
    %b: rf2: single rf pulse with 2 parts
    meas_name_suffix = '_MID0002_JB_RectPulse_a5';
    hc_name_suffix = '_MID0002_JB_RectPulse_a5';
    bc_name_suffix = '_MID0002_JB_RectPulse_a5_BC';
elseif subject_num == 3
    %c: rect rf pulses as reference
    meas_name_suffix = '_MID0003_JB_LIBRE2p2_a8_PERewinder';
    hc_name_suffix = '_MID0003_JB_LIBRE2p2_a8_PERewinder_HC';
    bc_name_suffix = '_MID0003_JB_LIBRE2p2_a8_PERewinder_BC';
elseif subject_num==4
    meas_name_suffix = '_MID0004_JB_LIBRE2p2_a8_woPERewinder';
    hc_name_suffix = '_MID0004_JB_LIBRE2p2_a8_woPERewinder_HC_state_clf';
    bc_name_suffix = '_MID0004_JB_LIBRE2p2_a8_woPERewinder_BC';
elseif subject_num==5
    meas_name_suffix = '_MID0005_JB_RectPulse_a5_PERewinder';
    hc_name_suffix = '_MID0005_JB_RectPulse_a5_PERewinder';
    bc_name_suffix = '_MID0005_JB_RectPulse_a5_PERewinder_BC';
elseif subject_num==6
    meas_name_suffix = '_MID0006_JB_RectPulse_a5_woPERewinder';
    hc_name_suffix = '_MID0006_JB_RectPulse_a5_woPERewinder';
    bc_name_suffix = '_MID0006_JB_RectPulse_a5_woPERewinder_BC';
else
    meas_name_suffix = '_MID0007_b_rfsp_adc';
    hc_name_suffix = '_MID0007_b_HC_rfsp_adc';
    bc_name_suffix = '_MID0007_b_BC_rfsp_adc';
end

meas_name = ['meas', meas_name_suffix];
hc_name = ['meas', hc_name_suffix];
bc_name = ['meas', bc_name_suffix];

measureFile = [datasetDir, meas_name,'.dat'];
bodyCoilFile = [datasetDir, bc_name,'.dat'];
arrayCoilFile = [datasetDir, hc_name,'.dat'];


%% Load and Configure Data
% Read data using the library's `createRawDataReader` function
% This readers makes the usage of Siemens and ISMRMRD files equivalent for
% the library
bodyCoilreader = createRawDataReader(bodyCoilFile, true);
% myTwix = mapVBVD_JH_for_monalisa(bodyCoilFile);
% bodyCoilreader.acquisitionParams.nSeg = 22;
% bodyCoilreader.acquisitionParams.nShot = 419;
% bodyCoilreader.acquisitionParams.nShot_off = 14;

arrayCoilReader = createRawDataReader(arrayCoilFile, true);
% arrayCoilReader.acquisitionParams.nSeg = 22;
% arrayCoilReader.acquisitionParams.nShot = 419;
% arrayCoilReader.acquisitionParams.nShot_off = 14;

% Ensure consistency in number o1f shot-off points
nShotOff = arrayCoilReader.acquisitionParams.nShot_off;



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
% Refine the sensitivity estimate with optimization
nIter = 5;
[C, x] = bmCoilSense_nonCart_secondary(y_surface, C_array_prime, y_ref, ...
                                       C_ref, Gn, Gu, Gut, ve, nIter, false);

% Display Results
bmImage(C);
%% Save C into the folder

saveCDirList = {strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_Binning/C/'),
    strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/C/')};


for idx = 2
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




