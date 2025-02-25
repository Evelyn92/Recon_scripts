%% Set up
clear;clc;


%% SPM to estimate the ref volume

addpath('/Users/cag/Documents/forclone/spm');

%% Start calculations

% 	1.	Load the 71 NIfTI volumes into SPM.
% 	2.	Perform realignment using spm_realign, which estimates motion parameters.
% 	3.	Apply the estimated motion correction using spm_reslice (optional).
% 	4.	Extract translation and rotation parameters from the rp_*.txt file.
% The reference volume is the mean of all volumes.

% Specify the directory containing the 71 NIfTI volumes
dataDir = '/Users/cag/Documents/Dataset/recon_results/250127_recon/Sub001/T1_LIBRE_Binning/output/mask_5s/nii';  

% Select all NIfTI files
niiFiles = spm_select('FPList', dataDir, '^.*\.nii$'); 

% Convert to cell array
files = cellstr(niiFiles);  

%% Initialize SPM
spm('defaults', 'FMRI');  
spm_jobman('initcfg');

%% Step 1: Estimate motion parameters (realignment)
matlabbatch = {};  % Initialize batch
matlabbatch{1}.spm.spatial.realign.estwrite.data = {files};  % Input files

% Set realignment parameters

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
%Defines the precision of motion estimation.

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.sep = 2;
%The distance (mm) between sampled points used for motion estimation

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.fwhm = 2;
% Applies a Gaussian filter (Full-Width Half-Maximum, FWHM) 
% to smooth the images before motion estimation. default 5mm


matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.rtm = 1; % 
	% •	0 (default): Each volume is realigned to the first volume.
	% •	1: All volumes are realigned to the mean image.

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.interp = 4;
% 	•	Defines how SPM interpolates voxel values when reslicing images.
% 	•	1: Nearest neighbor (fastest, but lowest quality).
% 	•	2 (default): Trilinear interpolation (good balance between quality and speed).
% 	•	4, 5: Higher-order spline interpolation (better quality but slower).

matlabbatch{1}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
	% •	Determines whether images are “wrapped” along each dimension.
	% •	[0 0 0]: No wrapping (default, used for most fMRI and structural images).
	% •	[1 1 0]: Wraps in X and Y directions (useful for phase-wrapped EPI images).

% Step 2: Apply motion correction (reslicing)
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.which = [2 1]; 
% Save realigned + mean image
%[2 1] → Saves both realigned images (with r* prefix) and the mean image (with mean* prefix).
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{1}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];

%% Run the batch
spm_jobman('run', matlabbatch);

%% Step 3: Extract motion parameters
motionParamsFile = fullfile(dataDir, 'rp_*.txt');  % Motion parameters file

motionParams = load(spm_select('FPList', dataDir, '^rp_.*\.txt$'));  % Load it

%% Step 4: Separate translation and rotation
translation = motionParams(:, 1:3);  % X, Y, Z translations (mm)
rotation = motionParams(:, 4:6);  % Pitch, Roll, Yaw (radians)

%
figure; 
subplot(2,1,1);
plot(motionParams(:, 1:3), 'LineWidth',2)
xlabel('Time (Volume)'); ylabel('translation (mm/)');
legend('Tx', 'Ty', 'Tz');
title('Translation Parameters');
grid on;

subplot(2,1,2);
plot(motionParams(:, 4:6), 'LineWidth',2)
xlabel('Time (Volume)'); ylabel('rotation (rad)');
legend( 'Rx', 'Ry', 'Rz');
title('Rotation Parameters');
grid on;
%% folder management
rfiles = dir(fullfile(dataDir, 'r*.nii')); % Find realigned NIfTI files
rdataDir = '/Users/cag/Documents/Dataset/recon_results/250127_recon/Sub001/T1_LIBRE_Binning/output/mask_5s/r_nii';
if ~isfolder(rdataDir)
    mkdir(rdataDir);
    disp(['The dir is created: ', rdataDir])
else
    disp(['The dir exist: ', rdataDir])
end
%
for i = 1:length(rfiles)
    oldName = fullfile(dataDir, rfiles(i).name);
    newName = fullfile(rdataDir,  rfiles(i).name); % Change prefix
    movefile(oldName, newName);
end

%% Visualization
idx = 69; % the slice index
o_img = files{idx};
r_img = fullfile(rdataDir, rfiles(idx).name);
mean_img = fullfile(dataDir, 'meanvolume_01.nii');
spm_check_registration(o_img, r_img);
% spm_check_registration(r_img, mean_img);
% spm_image('Display', r_img); 


      