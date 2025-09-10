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
addpath(genpath('/Users/cag/Documents/forclone/pulseq_v15'));
addpath(genpath('/Users/cag/Documents/forclone/monalisa'));



%% Initialize the directories and acquire the Coil
subject_num = 2;

%meas_MID00343_FID55397_yj_seq1.dat
% meas_MID00344_FID55398_yj_seq2_gdsp.dat
% 
% 

datasetDir = '/Users/cag/Documents/Dataset/datasets/250905/';
reconDir = '/Users/cag/Documents/Dataset/recon_results/250905/';
seqFolder = "/Users/cag/Documents/Dataset/datasets/250905/";
mask_note_list={'swap1_FA4_RF2_Shot2055','swap1_FA4_RF2_gdsp_Shot2055'};

mask_note = mask_note_list{subject_num};

if subject_num == 1
    meas_name_suffix = '_MID00343_FID55397_yj_seq1';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=2055;
    seqName = "yj_seq1_t1w_libre_main_TR6.2ms_TE3.6ms_swap1_FA4_RF2.seq";
elseif subject_num == 2
    meas_name_suffix = '_MID00344_FID55398_yj_seq2_gdsp';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=2055;
    seqName = "yj_seq1_t1w_libre_main_TR8.0ms_TE3.6ms_swap1_FA4_RF2_gdsp.seq";
end

meas_name = ['meas', meas_name_suffix];
hc_name = ['meas', hc_name_suffix];
bc_name = ['meas', bc_name_suffix];

measureFile = [datasetDir, meas_name,'.dat'];
bodyCoilFile = [datasetDir, bc_name,'.dat'];
arrayCoilFile = [datasetDir, hc_name,'.dat'];


% Load and Configure Data

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


%% ==============================================
% Warning: due to the memory limit, make sure the matrix size <=240
matrix_size = 240;  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./240;

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
disp(x0Path);
%

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


