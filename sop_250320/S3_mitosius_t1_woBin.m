clear; clc;
% =====================================================
% Author: Yiwei Jia
% Date: March 20
% ------------------------------------------------
% Coil sensitivity -> binning mask eMask -> [Mitosius]
% Update: the eyeGenerateBinning is replaced with eyeGenerateBinningWin
% Add a new txt file saved along side with the mask to log the details
% =====================================================

% change woBinning/Binning, change th75->th0, change eyeMask name!
%%
for subject_num = 5:6
subject_num_C = 3;
datasetDir = '/home/debi/yiwei/mreye_dataset/debug_0320/data/';
reconDir = '/home/debi/yiwei/recon_results/250320/';
% 
% ETDir = ['/home/debi/yiwei/recon_results/250314_seq_rfsp/', 'masks/'];
mask_note_list={'..','..', ...
    'JB_LIBRE2p2_a8_PERewinder', 'JB_LIBRE2p2_a8_woPERewinder', ...
    'JB_RectPulse_a5_PERewinder', 'JB_RectPulse_a5_woPERewinder'};
mask_note = mask_note_list{subject_num};


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




otherDir = [reconDir, strcat('/Sub00',num2str(subject_num),'/T1_LIBRE_woBinning/other/')];



% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end


saveCDirList = {strcat('/Sub00',num2str(subject_num_C),'/T1_LIBRE_Binning/C/'),
    strcat('/Sub00',num2str(subject_num_C),'/T1_LIBRE_woBinning/C/')};

%% Step 1: Load the Raw Data
autoFlag = true;  % Disable validation UI
reader = createRawDataReader(measureFile, autoFlag);
p = reader.acquisitionParams;
p.traj_type = 'full_radial3_phylotaxis';  % Trajectory type
% Acquisition from Bern need to change the following part!!
p.nShot_off = 14; % in case no validation UI
p.nShot = 2055; % in case no validation UI
p.nSeg = 22; % in case no validation UI
%%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(p);                       % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

%% Step 2: Load Coil Sensitivity Maps
% Load the coil sensitivity previously measured
saveCDir     = [reconDir,saveCDirList{2}];
CfileName = 'C.mat';
CfilePath = fullfile(saveCDir, CfileName);
load(CfilePath, 'C');  % Load sensitivity maps
disp(['C is loaded from:', CfilePath]);
% Adjust grid size for coil sensitivity maps
FoV = p.FoV;  % Field of View

% ==============================================
% Warning: due to the memory limit, all the voxel_size set on debi
% is always >= 1 to make sure the matrix size <=240
voxel_size = 2;
% So the mitosius saved on debi
% is the smaller than the full resolution.
% ===============================================
matrix_size = round(FoV/voxel_size);  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;
%
C = bmImResize(C, [48, 48, 48], N_u);
%% Step 3: Normalize the Raw Data
if N_u >240
    normalization = false;
else 
    normalization = true;
end
if normalization
    x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, N_u, dK_u); 
    %
    bmImage(x_tot)
    %
    temp_im = getimage(gca);  
    bmImage(temp_im); 
    temp_roi = roipoly; 
    normalize_val = mean(temp_im(temp_roi(:))); 
    % The normalize_val is super small, it is 5e-10, very small
    % again 3e-9
    % The value of one complex point is like: -0.0396 - 0.1162i
    disp('normalize_val')
    disp(normalize_val)
    y_tot(1,1,123)
end
% only once !!!!
if real(y_tot)<1
    if normalization
        y_tot = y_tot/normalize_val; 
        y_tot(1,1,123)
    else
        y_tot = y_tot/(3e-9); 
        y_tot(1,1,123)
    end
end



%% Set the folder for mitosius saving

mDir = [reconDir, strcat('/Sub00', num2str(subject_num)),'/T1_LIBRE_woBinning/mitosius/mask_', mask_note, '/'];

%% Prepare eye mask

eMaskFilePath = [otherDir,'eMask_woBin'];

eyeMask = load(eMaskFilePath); 
fields = fieldnames(eyeMask);  % Get the field names
firstField = fields{1};  % Get the first field name
eyeMask = eyeMask.(firstField);  % Access the first field's value
disp(eMaskFilePath)
disp('is loaded!')
% Eleminate the first segment of all the spokes for accuracies

%
size_Mask = size(eyeMask);
nbins = size_Mask(1);
eyeMask = reshape(eyeMask, [nbins, p.nSeg, p.nShot]); 
eyeMask(:, 1, :) = []; 

eyeMask(:, :, 1:p.nShot_off) = []; 
eyeMask = bmPointReshape(eyeMask); 


%% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot, t_tot, eyeMask); 
y = bmPermuteToCol(y); 
ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)


end





















