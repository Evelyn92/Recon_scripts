% =====================================================
% Author: Yiwei Jia
% Date: April 02
% ------------------------------------------------
% [Coil sensitivity] -> binning mask eMask -> Mitosius
% Update: this script is derived from Demo script
% by Mauro in Monalisa version Feb.5
% The old script has issue when running mask generation
% With readers, the param setting is more organized
% =====================================================
% 
% meas_MID00544_FID09845_v15_sub1_pre_HC.dat
% meas_MID00551_FID09852_v15_sub1_pre_BC.dat
% meas_MID00552_FID09853_v15_sub5_pre_HC.dat
% meas_MID00553_FID09854_v15_sub5_pre_BC.dat
% meas_MID00555_FID09856_v15_sub8_pre_HC.dat
% meas_MID00556_FID09857_v15_sub8_pre_BC.dat
% meas_MID00559_FID09860_v14_sub3_pre_HC.dat
% meas_MID00560_FID09861_v14_sub3_pre_BC.dat
% meas_MID00561_FID09862_v14_sub2_pre_HC.dat
% meas_MID00562_FID09863_v14_sub2_pre_BC.dat
% meas_MID00563_FID09864_v14_sub4_pre_HC.dat
% meas_MID00564_FID09865_v14_sub4_pre_BC.dat
% ----------------------------------
%--------- Acquired at Bern---------
% meas_MID00019_FID295654_v15_sub9_HC.dat
% meas_MID00022_FID295657_v15_sub9_BC.dat
% meas_MID00023_FID295658_v15_sub9_debug.dat
% meas_MID00024_FID295659_v15_sub10_HC.dat
% meas_MID00025_FID295660_v15_sub10_BC.dat
% meas_MID00026_FID295661_v15_sub10_debug.dat
% meas_MID00027_FID295662_v15_sub11_HC.dat
% meas_MID00028_FID295663_v15_sub11_BC.dat
% meas_MID00029_FID295664_v15_sub11_debug.dat
% meas_MID00030_FID295665_v15_sub12_HC.dat
% meas_MID00031_FID295666_v15_sub12_BC.dat
% meas_MID00032_FID295667_v15_sub12_debug.dat
%---repeat sub001 and sub002 as the baseline---
% for convenience: subject 13 and subject 14
% meas_MID00033_FID295668_v15_sub1_HC_Bern.dat
% meas_MID00034_FID295669_v15_sub1_BC_Bern.dat
% meas_MID00035_FID295670_v15_sub2_HC_Bern.dat
% meas_MID00037_FID295672_v15_sub2_BC_Bern.dat

%----------------------------------
clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));
%% Initialize the directories and acquire the Coil
subject_num = 14;
%
datasetDir = '/home/debi/yiwei/mreye_dataset/debug_0402/';
reconDir = '/home/debi/yiwei/recon_results/250402/';


if subject_num == 1
    %a: rf1: 2 rf pulses
    meas_name_suffix = '_MID00544_FID09845_v15_sub1_pre_HC';
    hc_name_suffix = '_MID00544_FID09845_v15_sub1_pre_HC';
    bc_name_suffix = '_MID00551_FID09852_v15_sub1_pre_BC';
elseif subject_num == 2
    % % meas_MID00561_FID09862_v14_sub2_pre_HC.dat
% meas_MID00562_FID09863_v14_sub2_pre_BC.dat
    %b: rf2: single rf pulse with 2 parts
    meas_name_suffix = '_MID00561_FID09862_v14_sub2_pre_HC';
    hc_name_suffix = '_MID00561_FID09862_v14_sub2_pre_HC';
    bc_name_suffix = '_MID00562_FID09863_v14_sub2_pre_BC';
elseif subject_num == 3

    meas_name_suffix = '_MID00559_FID09860_v14_sub3_pre_HC';
    hc_name_suffix = '_MID00559_FID09860_v14_sub3_pre_HC';
    bc_name_suffix = '_MID00560_FID09861_v14_sub3_pre_BC';
elseif subject_num==4
    meas_name_suffix = '_MID00563_FID09864_v14_sub4_pre_HC';
    hc_name_suffix = '_MID00563_FID09864_v14_sub4_pre_HC';
    bc_name_suffix = '_MID00564_FID09865_v14_sub4_pre_BC';
elseif subject_num==5
    meas_name_suffix = '_MID00552_FID09853_v15_sub5_pre_HC';
    hc_name_suffix = '_MID00552_FID09853_v15_sub5_pre_HC';
    bc_name_suffix = '_MID00553_FID09854_v15_sub5_pre_BC';
elseif subject_num==6
    meas_name_suffix = '_MID0006_JB_RectPulse_a5_woPERewinder';
    hc_name_suffix = '_MID0006_JB_RectPulse_a5_woPERewinder';
    bc_name_suffix = '_MID0006_JB_RectPulse_a5_woPERewinder_BC';
elseif subject_num == 8
    % _MID00556_FID09857_v15_sub8_pre_BC
    meas_name_suffix = '_MID00555_FID09856_v15_sub8_pre_HC';
    hc_name_suffix = '_MID00555_FID09856_v15_sub8_pre_HC';
    bc_name_suffix = '_MID00556_FID09857_v15_sub8_pre_BC';
elseif subject_num == 9
    meas_name_suffix = '_MID00019_FID295654_v15_sub9_HC';
    hc_name_suffix = '_MID00019_FID295654_v15_sub9_HC';
    bc_name_suffix = '_MID00022_FID295657_v15_sub9_BC';
elseif subject_num == 10
    meas_name_suffix = '_MID00024_FID295659_v15_sub10_HC';
    hc_name_suffix = '_MID00024_FID295659_v15_sub10_HC';
    bc_name_suffix = '_MID00025_FID295660_v15_sub10_BC';
elseif subject_num == 11
    meas_name_suffix = '_MID00027_FID295662_v15_sub11_HC';
    hc_name_suffix = '_MID00027_FID295662_v15_sub11_HC';
    bc_name_suffix = '_MID00028_FID295663_v15_sub11_BC';
elseif subject_num == 12
    meas_name_suffix = '_MID00030_FID295665_v15_sub12_HC';
    hc_name_suffix = '_MID00030_FID295665_v15_sub12_HC';
    bc_name_suffix = '_MID00031_FID295666_v15_sub12_BC';

elseif subject_num == 13
    meas_name_suffix = '_MID00033_FID295668_v15_sub1_HC_Bern';
    hc_name_suffix = '_MID00033_FID295668_v15_sub1_HC_Bern';
    bc_name_suffix = '_MID00034_FID295669_v15_sub1_BC_Bern';
elseif subject_num == 14
    meas_name_suffix = '_MID00035_FID295670_v15_sub2_HC_Bern';
    hc_name_suffix = '_MID00035_FID295670_v15_sub2_HC_Bern';
    bc_name_suffix = '_MID00037_FID295672_v15_sub2_BC';
    % bc_name_suffix = '_MID00148_FID291939_sub002_Bern_BC_replace_from0314';
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
bodyCoilreader.acquisitionParams.nSeg = 22;
bodyCoilreader.acquisitionParams.nShot = 419;
bodyCoilreader.acquisitionParams.nShot_off = 14;

arrayCoilReader = createRawDataReader(arrayCoilFile, true);
arrayCoilReader.acquisitionParams.nSeg = 22;
arrayCoilReader.acquisitionParams.nShot = 419;
arrayCoilReader.acquisitionParams.nShot_off = 14;

% Ensure consistency in number o1f shot-off points
nShotOff = arrayCoilReader.acquisitionParams.nShot_off;



%% Parameters
dK_u = [1, 1, 1] ./ arrayCoilReader.acquisitionParams.FoV;   % Cartesian grid spacing
N_u = [48, 48, 48];             % Adjust this value as needed
% Compute Trajectory and Volume Elements
[y_body, t, ve] = bmCoilSense_nonCart_data(bodyCoilreader, N_u);
y_surface = bmCoilSense_nonCart_data(arrayCoilReader, N_u);

% Compute the gridding matrices (subscript is a reminder of the result)
% Gn is from uniform to Non-uniform
% Gu is from non-uniform to Uniform
% Gut is Gu transposed
[Gn, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);
% Create Mask
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




