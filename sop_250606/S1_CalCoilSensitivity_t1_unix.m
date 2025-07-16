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

%% Initialize the directories and acquire the Coil
subject_num = 2;

datasetDir = '/home/debi/yiwei/mreye_dataset/250606/';
reconDir = '/home/debi/yiwei/recon_results/250606/';

% meas_MID00114_FID302856_YJ_seq1.dat
% meas_MID00115_FID302857_YJ_seq2.dat
% meas_MID00116_FID302858_YJ_seq3.dat
% meas_MID00117_FID302859_YJ_seq4.dat
% meas_MID00110_FID302852_JB_RectPulse_a5_woPERewinder_1.dat
% meas_MID00118_FID302860_JB_RectPulse_a5_woPERewinder_2_adjTRTE.dat


mask_note_list={'seq1_wrap_rfsp_no_gzsp','seq2_wrap_rfsp_no_gzsp_teRingdown',...
    'seq3_wrap_rfsp_gzsp_teRingdown','seq4_baseline_0603', ...
    'JB_rect_unchanged', 'JB_rect_yj_adj'};

mask_note = mask_note_list{subject_num};

if subject_num == 1
    meas_name_suffix = '_MID00114_FID302856_YJ_seq1';
    hc_name_suffix = '..';
    bc_name_suffix = '..';
    nShot=1000;

elseif subject_num == 2
    meas_name_suffix = '_MID00115_FID302857_YJ_seq2';
    hc_name_suffix = '..';
    bc_name_suffix = '..';
    nShot=1000;


elseif subject_num == 3
    meas_name_suffix = '_MID00116_FID302858_YJ_seq3';
    hc_name_suffix = '..';
    bc_name_suffix = '..';
    nShot=1000;

elseif subject_num == 4
    meas_name_suffix = '_MID00117_FID302859_YJ_seq4';
    hc_name_suffix = '..';
    bc_name_suffix = '..';
    nShot=1000;

elseif subject_num == 5
    meas_name_suffix = '_MID00110_FID302852_JB_RectPulse_a5_woPERewinder_1';
    hc_name_suffix = '..';
    bc_name_suffix = '..';
    nShot=2055;
else
    meas_name_suffix = '_MID00118_FID302860_JB_RectPulse_a5_woPERewinder_2_adjTRTE';
    hc_name_suffix = '..';
    bc_name_suffix = '..';
    nShot=1000;
end


meas_name = ['meas', meas_name_suffix];
% hc_name = ['meas', hc_name_suffix];
% bc_name = ['meas', bc_name_suffix];

measureFile = [datasetDir, meas_name,'.dat'];

%%
twix_1 = mapVBVD_JH_for_monalisa(measureFile);
%%
twix_2 = mapVBVD_JH_for_monalisa(measureFile);
%%
reader = createRawDataReader(measureFile, true);
%%
reader.acquisitionParams.nSeg = 22;
reader.acquisitionParams.nShot = nShot;
reader.acquisitionParams.nShot_off = 14;
reader.acquisitionParams.traj_type = 'full_radial3_phylotaxis';

% Ensure consistency in number o1f shot-off points
nShotOff = reader.acquisitionParams.nShot_off;

p = reader.acquisitionParams;
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
if (subject_num == 5) || (subject_num == 6)
    voxel_size = 1;
else
    voxel_size = 2;
end
% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = round(FoV/voxel_size);  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;
%%
% C = bmImResize(C, [48, 48, 48], N_u);
% x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u); 
% %
% bmImage(x_tot)
% %
% xtotDir = [reconDir, '/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/output/mask_',mask_note,'/'];
% if ~isfolder(xtotDir)
%     % If it doesn't exist, create it
%     mkdir(xtotDir);
%     disp(['Directory created: ', xtotDir]);
% else
%     disp(['Directory already exists: ', xtotDir]);
% end
% xtotPath = fullfile(xtotDir, 'x_tot_C.mat');
% % Save the x0 to the .mat file
% save(xtotPath, 'x_tot', '-v7.3');
% disp('xtot has been saved here:')
% disp(xtotPath)
%% ------
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
%%
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
%% Root mean square across the channels
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

clear;
