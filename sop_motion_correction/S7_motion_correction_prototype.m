clc;clear;
%% Some Flags
FlagSaveMask = false;
FlagReadrpTable = false;
FlagRot = true;
FlagTrans = true;
%% Get Path to the raw data
% Prompt user to select Siemens Raw Data file
basefolder = '/Users/cag/Documents/Dataset/datasets/250127_acquisition';
addpath(genpath(basefolder));
reconDir = '/Users/cag/Documents/Dataset/recon_results/250127_recon/';
resultsDir = fullfile(reconDir, 'temp_masks_2win_trial');  % Results folder
%
[fileName, filePath] = uigetfile({'*.dat', '(*.dat) Siemens Rawdatafile'}, ...
    'Select a Subject Rawdatafile', 'MultiSelect', 'off',basefolder);

if FlagReadrpTable
%%-- Start calculations
% Load the translation and rotation txt
rotTraslFile = "/Users/cag/Documents/Dataset/recon_results/250127_recon/Sub001/T1_LIBRE_Binning/output/mask_5s/nii_ref_first/rp_volume_01.txt";

% check if the rotTrasFile exists
if exist(fullfile(rotTraslFile),'file')
    disp(['The transform matrix txt file: ', rotTraslFile])
else
    disp('Warning!!! The transform matrix txt file does not exist!!!')
end

rpTable = importdata(rotTraslFile);
% Convert tables in to vectors 
TransParam=[rpTable(:,1)';rpTable(:,2)';rpTable(:,3)']; % translational components
RotParam=[rpTable(:,4)';rpTable(:,5)';rpTable(:,6)']; % rotational components
clear rpTable
end
%% Read the raw data .dat
brainScanFile = fullfile(filePath, fileName);  % Brain scan data file
% Create a reader object for the brain scan
brainCoilReader = createRawDataReader(brainScanFile, false);
% Extract acquisition parameters from the raw data
acquisitionParams = brainCoilReader.acquisitionParams;
nSeg = acquisitionParams.nSeg;  % Number of segments
nShotOff = acquisitionParams.nShot_off;  % Number of off shots (non-steady state)
nMeasuresPerShot = acquisitionParams.nSeg;  % Measurements per shot
nExcludeMeasures = nShotOff * nMeasuresPerShot;  % Total number of measurements to exclude
nLines = acquisitionParams.nLine;  % Total number of radial lines
nCh = acquisitionParams.nCh;
nShot = acquisitionParams.nShot;
NCol = acquisitionParams.N;
%% Calculate masks
if FlagSaveMask
% Define cost time (this is typically a fixed value depending on the scanner type)
costTime = 2.5;  % Siemens-specific, don't change unless known

% Extract timestamp in milliseconds (time of each acquisition)
timeStamp = double(acquisitionParams.timestamp);
timeStamp = timeStamp - min(timeStamp);  % Normalize timestamps to start from 0
timestampMs = timeStamp * costTime;  % Convert to milliseconds

%%%%%%% Bin the image with the first two 5s of sequence -> generate only 2 bin ->
%%%%%%%%% Teva_recon
% Define temporal window size in seconds (5 seconds in this case)
temporalWindowSec = 5;  
% Convert temporal window to milliseconds
temporalWindowMs = temporalWindowSec * 1000;
% Exclude non-steady-state measurements by considering the first few off shots
% Adjust the start time to account for non-steady-state shots
startTime = timestampMs(nExcludeMeasures + 1);
% Calculate the total duration of valid data
nMasks = 2; %here I set nBin=2 just for testing
totalDuration = temporalWindowMs*nMasks;
endTime = startTime+totalDuration;

% Initialize the mask matrix with logical false
mask = false(nMasks, nLines);
% Fill the masks: Set to true for measurements that fall within the current temporal window
for i = 1:nMasks
    % Define the start and end of the current time window
    windowStart = startTime + (i - 1) * temporalWindowMs;
    windowEnd = windowStart + temporalWindowMs;
    
    % Create the mask for the current window (True for measurements within the window)
    singlemask = (timestampMs >= windowStart) & (timestampMs < windowEnd);
    
    % Exclude the non-steady-state lines (e.g., SI projection or other artifacts)
    for K = 0:floor(nLines / nMeasuresPerShot)
        idx = 1 + K * nSeg;
        if idx <= nLines
            singlemask(idx) = false;  % Set non-steady-state lines to false
        end
    end
    
    % Assign the mask to the binning matrix
    mask(i, :) = singlemask;
end
% Save Results
% Generate a timestamped filename for saving
if ~isfolder(resultsDir)
    % If it doesn't exist, create it
    mkdir(resultsDir);
    disp(['Directory created: ', resultsDir]);
end
saveName = fullfile(resultsDir, 'win2Binning.mat');

save(saveName, 'mask');
disp(['Sequential bins saved to: ', saveName]);

end
%% Test with observable translation and rotation
% FakeTrans = repmat([5; 2; 0], 1, 71);
FakeTrans = repmat([0; 48; 0], 1, 71);
FakeRot = repmat([0;0;0],1, 71);

%Here add minus here for inverse direction
TransParam=FakeTrans;
RotParam=FakeRot;
              
%% bmTraj for the full traj calculation
% In order to get the nLine of traj=total meas
% here I set selfNav_flag = False
% nShot_off = 0;
p = acquisitionParams;
p.traj_type = 'full_radial3_phylotaxis';  % Trajectory type
t_tot = bmTraj(p);                       % Compute filtered trajectory

% Note: any change in p will change the acquisitionParams
% need to put the right value back after computing the full trajectory.
acquisitionParams.nShot_off = nShotOff;


%% Rotation
% convert the Euler angle in radians to 3 by 3 rotation matrix  
if FlagRot
    iLine = 40;
    Rot_mat = euler_to_rotation(RotParam(:, iLine));
    
    t_tot_reshape =    reshape(t_tot, 3, []);    
    t_tot_rotated = Rot_mat * t_tot_reshape; 
    t_tot_rotated = reshape(t_tot_rotated, 3, NCol, []);

    t_tot_clean = t_tot_rotated;
    % newkx=reshape(t_tot_rotated(1,:,:),[param.Np,param.Nseg,param.Nshot]);
    % newky=reshape(t_tot_rotated(2,:,:),[param.Np,param.Nseg,param.Nshot]);
    % newkz=reshape(t_tot_rotated(3,:,:),[param.Np,param.Nseg,param.Nshot]);
else
    t_tot_clean =  t_tot;


end
%% Translation
if FlagTrans
    kx = t_tot(1,:,:);
    ky = t_tot(2,:,:);
    kz = t_tot(3,:,:);
    % Note: here we calculate with the original trajectory, not the new one from rotated...           
    
    tempx = kx.*TransParam(1,iLine);       
    tempy = ky.*TransParam(2,iLine);   
    tempz = kz.*TransParam(3,iLine);   
    PhaseOffset = exp(-2*pi*1i*(tempx + tempy + tempz));
    % add the PhaseOffset to y
    y_tot = brainCoilReader.readRawData(true, true); % DO filter nshotoff and SI
    y_tot_phase = y_tot .* PhaseOffset;
    y_tot_clean = y_tot_phase;

               
else
    y_tot_clean = brainCoilReader.readRawData(true, true); % DO NOT filter nshotoff and SI

end

%% mitosius

correction_mitosius(t_tot_clean, y_tot_clean, reconDir, resultsDir, acquisitionParams);

% reconstruction
x0 = correction_reconstruction(reconDir);

 %%
function R = euler_to_rotation(R_array)
    % Compute rotation matrices for X, Y, Z axes
    Rx= R_array(1); 
    Ry= R_array(2); 
    Rz= R_array(3);

    Rx_mat = [1  0       0      ;
              0  cos(Rx) -sin(Rx);
              0  sin(Rx)  cos(Rx)];

    Ry_mat = [cos(Ry)  0  sin(Ry);
              0        1  0      ;
             -sin(Ry)  0  cos(Ry)];

    Rz_mat = [cos(Rz) -sin(Rz)  0;
              sin(Rz)  cos(Rz)  0;
              0        0        1];

    % Combine rotations: R = Rz * Ry * Rx
    R = Rz_mat * Ry_mat * Rx_mat;
end


function correction_mitosius(t_tot_clean, y_tot_clean, reconDir, resultsDir, acquisitionParams)
% Prepare the input variable
t_tot = t_tot_clean;
y_tot = y_tot_clean;
p = acquisitionParams;
mask_note = 'mc_2bin_5s';
mDir = [reconDir, '/Sub001/T1_LIBRE_Binning/mitosius/mask_', mask_note, '/'];


% Load the binning mask
tempMaskFilePath = fullfile(resultsDir, 'win2Binning.mat');
tempMask = load(tempMaskFilePath); 
fields = fieldnames(tempMask);  % Get the field names
firstField = fields{1};  % Get the first field name
tempMask = tempMask.(firstField);  % Access the first field's value
disp(tempMaskFilePath)
disp('is loaded!')

% Prepare the binning mask
size_Mask = size(tempMask);
nbins = size_Mask(1);
tempMask = reshape(tempMask, [nbins, p.nSeg, p.nShot]); 
tempMask(:, 1, :) = []; 
p=acquisitionParams;
tempMask(:, :, 1:p.nShot_off) = []; 
tempMask = bmPointReshape(tempMask); 

% Load the coil sensitivity previously measured
saveCDir     = [reconDir,'/Sub001/T1_LIBRE_Binning/C/']; CfileName = 'C.mat';
CfilePath = fullfile(saveCDir, CfileName); load(CfilePath, 'C');  % Load sensitivity maps
disp(['C is loaded from:', CfilePath]);
% Adjust grid size for coil sensitivity maps
FoV = p.FoV;  % Field of View
voxel_size = 2;
matrix_size = round(FoV/voxel_size);  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;
C = bmImResize(C, [48, 48, 48], N_u);

% Normalize the Raw Data
if real(y_tot)<1
    y_tot = y_tot/(3e-9); 
    y_tot(1,1,123)
    disp('y_tot is normalized by 3e-9.')
else
    disp("it is already normalized!")
end

[y, t] = bmMitosis(y_tot, t_tot, tempMask); 
y = bmPermuteToCol(y); 
ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)

end


function x0 = correction_reconstruction(reconDir)
mask_note = 'mc_2bin_5s';
mDir = [reconDir, '/Sub001/T1_LIBRE_Binning/mitosius/mask_', mask_note, '/'];
CfileName = 'C.mat'; saveCDir = [reconDir, '/Sub001/T1_LIBRE_Binning/C/'];
CfilePath = fullfile(saveCDir, CfileName); load(CfilePath); 

y   = bmMitosius_load(mDir, 'y'); 
t   = bmMitosius_load(mDir, 't'); 
ve  = bmMitosius_load(mDir, 've');
disp('Mitosius has been loaded!')

Matrix_size = 120;
ReconFov = 240; %mm
N_u     = [Matrix_size, Matrix_size, Matrix_size]; % Matrix size: Size of the Virtual cartesian grid in the fourier space (regridding)
n_u     = [Matrix_size, Matrix_size, Matrix_size]; % Image size (output)
dK_u    = [1, 1, 1]./ReconFov; % Spacing of the virtual cartesian grid
nFr     = size(y,1); 

C = bmImResize(C, [48, 48, 48], N_u);

x0 = cell(nFr, 1);
for i = 1:nFr
    x0{i} = bmMathilda(y{i}, t{i}, ve{i}, C, N_u, n_u, dK_u, [], [], [], []);
end

%
bmImage(x0);


end