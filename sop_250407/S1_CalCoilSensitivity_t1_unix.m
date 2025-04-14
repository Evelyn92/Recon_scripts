% =====================================================
% Author: Yiwei Jia
% Date: April 07
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
subject_num = 3;
%
datasetDir = '/home/debi/yiwei/mreye_dataset/250407/';
reconDir = '/home/debi/yiwei/recon_results/250407/';

if subject_num == 1
    meas_name_suffix = '_MID00333_FID11707_sub1_debug_brain';
    hc_name_suffix = '_MID00334_FID11708_sub1_HC_brain';
    bc_name_suffix = '_MID00341_FID11715_sub1_BC_brain';

elseif subject_num == 2
    meas_name_suffix = '_MID00342_FID11716_sub2_debug_brain';
    hc_name_suffix = '_MID00343_FID11717_sub2_HC_brain';
    bc_name_suffix = '_MID00344_FID11718_sub2_BC_brain';

elseif subject_num == 3
    meas_name_suffix = '_MID00345_FID11719_sub3_debug_brain';
    hc_name_suffix = '_MID00346_FID11720_sub3_HC_brain';
    bc_name_suffix = '_MID00347_FID11721_sub3_BC_brain';
elseif subject_num == 4
    meas_name_suffix = '_MID00309_FID11683_sub1_HC';
    hc_name_suffix = '_MID00309_FID11683_sub1_HC';
    bc_name_suffix = '_MID00316_FID11690_sub1_BC';

elseif subject_num == 5
    meas_name_suffix = '_MID00317_FID11691_sub2_HC';
    hc_name_suffix = '_MID00317_FID11691_sub2_HC';
    bc_name_suffix = '_MID00318_FID11692_sub2_BC';
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
bodyCoilreader.acquisitionParams.traj_type = 't5';

arrayCoilReader = createRawDataReader(arrayCoilFile, true);
arrayCoilReader.acquisitionParams.nSeg = 22;
arrayCoilReader.acquisitionParams.nShot = 419;
arrayCoilReader.acquisitionParams.nShot_off = 14;
arrayCoilReader.acquisitionParams.traj_type = 't5';

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
[C1, x] = bmCoilSense_nonCart_secondary(y_surface, C_array_prime, y_ref, ...
                                       C_ref, Gn, Gu, Gut, ve, nIter, false);

% Display Results
bmImage(C1);
%% 
C = C1;
for iCh = 1:size(C1,4)
    C(:,:,:,iCh) = flip(permute(C1(:,:,:,iCh), [2 1 3]),2);
    C(:,:,:,iCh) = flip(permute(C(:,:,:,iCh), [2 1 3]),2);
end
bmImage(C)
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




%%

p = arrayCoilReader.acquisitionParams;
p.traj_type = 'rot_x_full_radial3_phylotaxis';  % Trajectory type
% Acquisition from Bern need to change the following part!!
p.nShot_off = 14; % in case no validation UI
p.nShot = nShot; % in case no validation UI
p.nSeg = 22; % in case no validation UI
%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(p);                       % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

% Adjust grid size for coil sensitivity maps
FoV = p.FoV;  % Field of View

% ==============================================
% Warning: due to the memory limit, all the voxel_size set on debi
% is always >= 1 to make sure the matrix size <=240
voxel_size = 4;
% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = round(FoV/voxel_size);  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u); 