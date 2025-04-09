clear; clc;
% =====================================================
% Author: Yiwei Jia
% Date: April 03
% ------------------------------------------------
% Coil sensitivity -> binning mask eMask -> [Mitosius]
% Update: the eyeGenerateBinning is replaced with eyeGenerateBinningWin
% Add a new txt file saved along side with the mask to log the details
% =====================================================

% change woBinning/Binning, change th75->th0, change eyeMask name!
%%
subject_array = 14;
% for subject_num = subject_array
subject_num = subject_array;
subject_num_C = subject_array;

datasetDir = '/home/debi/yiwei/mreye_dataset/debug_0402/';
reconDir = '/home/debi/yiwei/recon_results/250402/';
% 
% ETDir = ['/home/debi/yiwei/recon_results/250314_seq_rfsp/', 'masks/'];
mask_note_list={'libre_wo_rfsp','rect_wo_rfsp', ...
    'libre_rfsp_qua50_adc_ph', 'rect_rfsp_qua50_adc_ph', ...
    '..','..', ...
    '..','..', ...
    'libre_rfsp_qua50_adc_ph','libre_rfsp_qua-50_adc_-ph', ...
    'libre_rfsp_qua50_adc_-ph','libre_rfsp_qua50_adc_0', ...
    'libre_wo_rfsp', 'rect_wo_rfsp'};

mask_note = mask_note_list{subject_num};

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
    % bc_name_suffix = '_MID00037_FID295672_v15_sub2_BC_Bern';
    bc_name_suffix = '_MID00148_FID291939_sub002_Bern_BC_replace_from0314';
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
p.nShot = 419; % in case no validation UI
p.nSeg = 22; % in case no validation UI
%
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
voxel_size = 4;
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
    x0=x_tot;
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


% end





















