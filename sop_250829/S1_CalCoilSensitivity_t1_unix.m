% =====================================================
% Author: Yiwei Jia
% Date: June 05
% ------------------------------------------------
% [Coil sensitivity] -> binning mask eMask -> Mitosius
% Update: this script is derived from Demo script
% by Mauro in Monalisa version Feb.5
% The old script has issue when running mask generation
% With readers, the param setting is more organized
% =====================================================

clc;
addpath(genpath('/Users/cag/Documents/forclone/Recon_scripts'));
addpath(genpath('/Users/cag/Documents/forclone/pulseq'));
addpath(genpath('/Users/cag/Documents/forclone/monalisa'));
% meas_MID00157_FID314156_JB_LIBRE2p2_a8_woPERewinder.dat
% meas_MID00158_FID314157_t1w_swap_0.dat
% meas_MID00159_FID314158_t1w_swap_1.dat
% meas_MID00160_FID314159_JB_t2w_LIBRE2p2_a8_woPERewinder.dat
% meas_MID00184_FID314183_t2_22_hs.dat

% meas_MID00176_FID314175_JB_LIBRE2p2_t2w.dat small size
% meas_MID00177_FID314176_t2_22_hs.dat small


%% Initialize the directories and acquire the Coil
subject_num = 6;
use_C = 0;
%
datasetDir = '/Users/cag/Documents/Dataset/datasets/250829/';
reconDir = '/Users/cag/Documents/Dataset/recon_results/250829/';

mask_note_list={'swap1_FA4_RF2','swap1_FA8_RF2', ...
    'swap1_FA4_RF2_freqPos', 'swap1_FA16_RF2.seq', 'swap0_FA4_RF2.seq', 'idea_240'};

mask_note = mask_note_list{subject_num};

if subject_num == 1
    meas_name_suffix = '_MID00205_FID315299_t1w_seq1';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;
elseif subject_num == 2
    meas_name_suffix = '_MID00206_FID315300_t1w_seq2';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;

elseif subject_num == 3
    meas_name_suffix = '_MID00211_FID315305_t1w_seq3';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;

elseif subject_num == 4
    meas_name_suffix = '_MID00210_FID315304_t1w_seq4';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;
elseif subject_num == 5
    meas_name_suffix = '_MID00207_FID315301_t1w_seq5';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;

elseif subject_num == 6
    meas_name_suffix = '_MID00209_FID315303_JB_LIBRE2p2_a8_woPERewinder';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;
end


meas_name = ['meas', meas_name_suffix];
hc_name = ['meas', hc_name_suffix];
bc_name = ['meas', bc_name_suffix];

measureFile = [datasetDir, meas_name,'.dat'];
bodyCoilFile = [datasetDir, bc_name,'.dat'];
arrayCoilFile = [datasetDir, hc_name,'.dat'];


%% Load and Configure Data
if use_C
% Read data using the library's `createRawDataReader` function
% This readers makes the usage of Siemens and ISMRMRD files equivalent for
% the library
bodyCoilreader = createRawDataReader(bodyCoilFile, true);
% myTwix = mapVBVD_JH_for_monalisa(bodyCoilFile);
bodyCoilreader.acquisitionParams.nSeg = 22;
bodyCoilreader.acquisitionParams.nShot = 419;
bodyCoilreader.acquisitionParams.nShot_off = 14;
bodyCoilreader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
%
arrayCoilReader = createRawDataReader(arrayCoilFile, true);
arrayCoilReader.acquisitionParams.nSeg = 22;
arrayCoilReader.acquisitionParams.nShot = 419;
arrayCoilReader.acquisitionParams.nShot_off = 14;
arrayCoilReader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';

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
% for iCh = 1:size(C1,4)
%     C(:,:,:,iCh) = flip(permute(C1(:,:,:,iCh), [2 1 3]),2);
%     C(:,:,:,iCh) = flip(permute(C(:,:,:,iCh), [2 1 3]),2);
% end
bmImage(C)
% Save C into the folder

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


else

%
reader = createRawDataReader(measureFile, false);
reader.acquisitionParams.nSeg = 22;

reader.acquisitionParams.nShot_off = 14;
% reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';


% Ensure consistency in number o1f shot-off points
nShotOff = reader.acquisitionParams.nShot_off;

p = reader.acquisitionParams;

% Acquisition from Bern need to change the following part!!
p.nShot_off = 14; % in case no validation UI
p.nShot = nShot; % in case no validation UI
p.nSeg = 22; % in case no validation UI

%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(p);                       % Compute trajectory

ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./240;
end

if use_C
%
C = bmImResize(C, [48, 48, 48], N_u);
x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u); 
%
bmImage(x_tot)
%
xtotDir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/output/mask_',mask_note,'/'];
if ~isfolder(xtotDir)
    % If it doesn't exist, create it
    mkdir(xtotDir);
    disp(['Directory created: ', xtotDir]);
else
    disp(['Directory already exists: ', xtotDir]);
end
xtotPath = fullfile(xtotDir, 'x_tot_C.mat');
% Save the x0 to the .mat file
save(xtotPath, 'x_tot', '-v7.3');
disp('xtot has been saved here:')
disp(xtotPath)

else
% ------
nCh = size(y_tot, 1);
nCh
nFr = 1;
x0 = cell(nCh, 1);
for i = 1:nFr
    for iCh = 1:nCh
    x0{iCh} = bmMathilda(y_tot(iCh,:), t_tot, ve_tot, [], N_u, N_u, dK_u, [], [], [], []);
    disp(['Processing channel: ', num2str(iCh),'/', num2str(nCh)])
   
    end
end

%
bmImage(x0);

%
x0Dir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/output/mask_',mask_note,'/'];
 
if ~isfolder(x0Dir)
    % If it doesn't exist, create it
    mkdir(x0Dir);
    disp(['Directory created: ', x0Dir]);
else
    disp(['Directory already exists: ', x0Dir]);
end
x0Path = fullfile(x0Dir, 'x0.mat');
% Save the x0 to the .mat file
save(x0Path, 'x0', '-v7.3');
disp('x0 has been saved here:')
disp(x0Path)
% Root mean square across the channels
% Initialize an array to store sum of squared images
[nx, ny, nz] = size(x0{1});  % Get the dimensions (240,240,240)
numCoils = numel(x0);  % Number of coils (20)

sum_of_squares = zeros(nx, ny, nz, 'single');  % Preallocate in single precision

% Compute sum of squared images
for coil = 1:numCoils
    sum_of_squares = sum_of_squares + abs(x0{coil}).^2;
end

% Compute the root mean square (RMS)
xrms = sqrt(sum_of_squares / numCoils);  % Normalize by the number of coils
xrmsPath = fullfile(x0Dir, 'xrms.mat');

% Save the x0 to the .mat file
save(xrmsPath, 'xrms', '-v7.3');
disp('xrmsPath has been saved here:')
disp(xrmsPath)
%%
bmImage(xrms)
end

