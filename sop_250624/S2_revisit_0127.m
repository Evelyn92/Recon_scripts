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

clear; clc;
addpath(genpath('/home/debi/yiwei/forclone/Recon_scripts'));
addpath(genpath('/home/debi/yiwei/forclone/pulseq'));
% meas_MID00387_FID34124_pulseqv15.dat     
% meas_MID00388_FID34125_pulseqv15.dat    
% meas_MID00389_FID34126_pulseqv15.dat
% meas_MID00390_FID34127_pulseqv15.dat
% meas_MID00401_FID34138_pulseqv15_HC_BC.dat 
% meas_MID00400_FID34137_pulseqv15_BC_BC.dat 



%% Initialize the directories and acquire the Coil
subject_num = 5;

%
datasetDir = '/home/debi/yiwei/mreye_dataset/250624/';
reconDir = '/home/debi/yiwei/recon_results/250624/';

mask_note_list={'seq1_t1w_gre_TR6p2','seq1_t1w_libre_TR6p2', ...
    'seq2_t1w_gre_TR8p01','seq2_t1w_libre_TR8p01','revisit_0127'

    };

mask_note = mask_note_list{subject_num};

if subject_num == 1
    meas_name_suffix = '_MID00387_FID34124_pulseqv15';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;

elseif subject_num == 2
    meas_name_suffix = '_MID00388_FID34125_pulseqv15';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;

elseif subject_num == 3
    meas_name_suffix = '_MID00389_FID34126_pulseqv15';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;
elseif subject_num == 4
    meas_name_suffix = '_MID00390_FID34127_pulseqv15';
    hc_name_suffix = ' ';
    bc_name_suffix = ' ';
    nShot=1000;

elseif subject_num == 5
    datasetDir = '/home/debi/yiwei/mreye_dataset/250127_acquisition/';
    meas_name_suffix = '_MID00332_FID214628_BEAT_LIBREon_eye_(23_09_24)_sc_trigger';
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

%
reader = createRawDataReader(measureFile, true);
if subject_num == 5
    reader.acquisitionParams.nShot_off = 14;
else
    reader.acquisitionParams.nSeg = 22;
    reader.acquisitionParams.nShot = 1000;
    reader.acquisitionParams.nShot_off = 14;
end




% Ensure consistency in number o1f shot-off points
nShotOff = reader.acquisitionParams.nShot_off;

p = reader.acquisitionParams;
if subject_num == 2
    p.traj_type = 'pulseq';
    p.pulseqTrajFile_name = "/home/debi/yiwei/mreye_dataset/250624/pulseq/seq1_t1w_libre_debugTR6p2_rfdelay_gain_0.seq";
else
    p.traj_type = 'full_radial3_phylotaxis';
end

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
voxel_size = 1;

% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = round(FoV/voxel_size);  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

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


