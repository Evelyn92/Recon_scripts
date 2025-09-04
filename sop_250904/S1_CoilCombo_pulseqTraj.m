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



%% Initialize the directories and acquire the Coil
subject_num = 1;

%
datasetDir = '/Users/cag/Documents/Dataset/datasets/250829/';
reconDir = '/Users/cag/Documents/Dataset/recon_results/250829/';
seqFolder = "/Users/cag/Documents/forclone/pulseq4mreye/dev/libre_3d_radial/output/0904_t1w/";
mask_note_list={'swap1_FA4_RF2_Shot1000','swap1_FA4_RF2_overlap0_Shot987', ...
    'swap1_FA4_RF2_overlap1_Shot987', 'swap1_FA4_RF2_overlap1_Shot987_BW2', 'swap1_FA4_RF2_overlap1_Shot987_FOV240_osR1'};

mask_note = mask_note_list{subject_num};

if subject_num == 1
    meas_name_suffix = '_MID00205_FID315299_t1w_seq1';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;
    seqName = "seq1_t1w_libre_part_TR6.2ms_TE3.6ms_swap1_FA4_RF2_Shot1000.seq";
elseif subject_num == 2
    meas_name_suffix = '_MID00206_FID315300_t1w_seq2';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=987;
    seqName = "seq2_t1w_libre_part_TR6.2ms_TE3.6ms_swap1_FA4_RF2_overlap0_Shot987.seq";
elseif subject_num == 3
    meas_name_suffix = '_MID00211_FID315305_t1w_seq3';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=987;
    seqName = "seq3_t1w_libre_part_TR6.2ms_TE3.6ms_swap1_FA4_RF2_overlap1_Shot987.seq";
elseif subject_num == 4
    meas_name_suffix = '_MID00210_FID315304_t1w_seq4';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=987;
    seqName = "seq4_t1w_libre_part_TR6.2ms_TE3.6ms_swap1_FA4_RF2_overlap1_Shot987_BW2.seq";
elseif subject_num == 5
    meas_name_suffix = '_MID00207_FID315301_t1w_seq5';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=987;
    seqName = "seq5_t1w_libre_part_TR6.2ms_TE3.6ms_swap1_FA4_RF2_overlap1_Shot987_FOV240_osR1.seq";
end


meas_name = ['meas', meas_name_suffix];
hc_name = ['meas', hc_name_suffix];
bc_name = ['meas', bc_name_suffix];

measureFile = [datasetDir, meas_name,'.dat'];
bodyCoilFile = [datasetDir, bc_name,'.dat'];
arrayCoilFile = [datasetDir, hc_name,'.dat'];


%% Load and Configure Data

reader = createRawDataReader(measureFile, false);
% Acquisition from Bern need to manually define the following part!!
reader.acquisitionParams.nSeg = 22;
reader.acquisitionParams.nShot = nShot; % in case no validation UI
reader.acquisitionParams.nShot_off = 14;
% reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';
reader.acquisitionParams.traj_type = 'pulseq';
reader.acquisitionParams.pulseqTrajFile_name = seqFolder + ...
    seqName;

% Ensure consistency in number o1f shot-off points
nShotOff = reader.acquisitionParams.nShot_off;
%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(reader.acquisitionParams);                       % Compute trajectory


ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

% Adjust grid size for coil sensitivity maps


% ==============================================
% Warning: due to the memory limit, make sure the matrix size <=240
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./p.FoV;

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

bmImage(xrms)


