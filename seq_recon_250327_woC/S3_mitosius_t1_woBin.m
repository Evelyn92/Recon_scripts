clear; clc;
% =====================================================
% Author: Yiwei Jia
% Date: March 27
% ------------------------------------------------
% Coil sensitivity -> binning mask eMask -> [Mitosius]
% Update: the eyeGenerateBinning is replaced with eyeGenerateBinningWin
% Add a new txt file saved along side with the mask to log the details
% =====================================================

% change woBinning/Binning, change th75->th0, change eyeMask name!
%%
subject_num = 2;
datasetDir = '/home/debi/yiwei/mreye_dataset/250324/v14_bern/';
reconDir = '/home/debi/yiwei/recon_results/250324_woC/';
mask_note_list = {'noC_libre','noC_rect'};
mask_note = mask_note_list{subject_num};

%%

if subject_num == 1
    bodyCoilFile     = [datasetDir, ' '];
    arrayCoilFile    = [datasetDir, ' '];
    measureFile = fullfile(datasetDir,'meas_MID00147_FID294264_libre_failed.dat');

elseif subject_num == 2
    bodyCoilFile     = [datasetDir, ' '];
    arrayCoilFile    = [datasetDir, ' '];
    measureFile = fullfile(datasetDir,'meas_MID00148_FID294265_rect_fine.dat');

elseif subject_num == 3
    bodyCoilFile = [datasetDir, ' '];
    arrayCoilFile = [datasetDir, ' '];
    measureFile = [datasetDir, ' '];

end


otherDir = [reconDir, '/Sub00', num2str(subject_num),'/T1_LIBRE_woBinning/other/'];

% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end


if subject_num == 1
    saveCDirList = {'/Sub001/T1_LIBRE_Binning/C/','/Sub001/T1_LIBRE_woBinning/C/'};
elseif subject_num == 2
    saveCDirList = {'/Sub002/T1_LIBRE_Binning/C/','/Sub002/T1_LIBRE_woBinning/C/'};
elseif subject_num == 3
    saveCDirList = {'/Sub003/T1_LIBRE_Binning/C/','/Sub003/T1_LIBRE_woBinning/C/'};
else
    saveCDirList = {'/Sub004/T1_LIBRE_Binning/C/','/Sub004/T1_LIBRE_woBinning/C/'};
end

%% Step 1: Load the Raw Data
autoFlag = false;  % Disable validation UI
reader = createRawDataReader(measureFile, autoFlag);
p = reader.acquisitionParams;
p.traj_type = 'full_radial3_phylotaxis';  % Trajectory type
% Acquisition from Bern need to change the following part!!
p.nShot_off = 5; % in case no validation UI
p.nShot = 419; % in case no validation UI
p.nSeg = 22; % in case no validation UI
%%
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(p);                       % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements

%% Step 2: Load Coil Sensitivity Maps


FoV = p.FoV/2; % due to setting from pulseq..

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
C = [];
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
% plot3(squeeze(t_tot(1,end,1:40)),squeeze(t_tot(2,end,1:40)), squeeze(t_tot(3,end,1:40)))
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


if subject_num == 1
    mDir = [reconDir, '/Sub001/T1_LIBRE_woBinning/mitosius/mask_', mask_note, '/'];
elseif subject_num == 2
    mDir = [reconDir,'/Sub002/T1_LIBRE_woBinning/mitosius/mask_', mask_note, '/'];
elseif subject_num == 3
    mDir = [reconDir, '/Sub003/T1_LIBRE_woBinning/mitosius/mask_', mask_note, '/'];
else
    mDir = [reconDir, '...'];
end

%%
mask = ones(1,p.nLine);
mask = mask>0;
mask = reshape(mask, [1, p.nSeg, p.nShot]); 
mask(:, 1, :) = []; 

mask(:, :, 1:p.nShot_off) = []; 
mask = bmPointReshape(mask); 




%% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot, t_tot, mask); 
y = bmPermuteToCol(y); 
ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)






















