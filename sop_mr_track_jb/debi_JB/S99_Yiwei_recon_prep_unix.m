clc;
% addpath(genpath("/media/sinf/1,0 TB Disk/Backup/Recon_fork"));
addpath(genpath('/media/sinf/1,0 TB Disk/Backup/Recon_fork'));
%%
% Coil sensitivity -> binning mask cMask -> Mitosius
% ------------------------------------------------------------------
% ------------------------------------------------------------------

% ------------------------------------------------------------------
% ------------------------------------------------------------------
%% Initialize the directories and acquire the Coil
subject_num = 4;
%
reconDir = '/media/sinf/1,0 TB Disk/240922_480/';
if subject_num == 1
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/RawData';
    ETDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/ET_EDF';
elseif subject_num == 2
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/MR_EYE_Subj02/RawData';
    ETDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/ET_EDF';
elseif subject_num == 3
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub003/230928_anatomical_MREYE_study/MR_EYE_Subj03/RawData';
    ETDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub003/230928_anatomical_MREYE_study/ET_EDF';
else
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub004/230923_anatomical_MREYE_study/MR_EYE_Subj04/RawData';
    ETDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub004/230923_anatomical_MREYE_study/ET_EDF';
end

if subject_num == 1
    bodyCoilFile     = [datasetDir, '/meas_MID00469_FID57935_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00470_FID57936_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00453_FID57919_BEAT_LIBREon_eye.dat'];
elseif subject_num == 2
    bodyCoilFile     = [datasetDir, '/meas_MID00357_FID56836_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00358_FID56837_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00342_FID56821_BEAT_LIBREon_eye.dat'];
elseif subject_num == 3
    bodyCoilFile = [datasetDir, '/meas_MID00312_FID57778_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile = [datasetDir, '/meas_MID00313_FID57779_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00297_FID57763_BEAT_LIBREon_eye.dat'];
else
    bodyCoilFile = [datasetDir, '/meas_MID00349_FID57815_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile = [datasetDir, '/meas_MID00350_FID57816_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00333_FID57799_BEAT_LIBREon_eye.dat'];
end
%% Start coil estimation
disp('-------------------------------------')
disp('Start Coil Estimation')
disp('-------------------------------------')
% Read metadata from the twix
bmTwix_info(bodyCoilFile)
% bmTwix_info(arrayCoilFile)
%%
inspect_raw_data = true;
if inspect_raw_data
    bmTwix_info(measureFile);
end

%% The input initialization is not automatizable, since the name/content of
% the variables depends on the sequence programmed. For the moment we just
% read and print the twixinfo and you need to adjut the parameters yourself

% Maybe: write a function for automate ISMR rawData format(Standard) reading.
% All trajectory information, to generate the trajectory. 
N            = 480; 
nSeg         = 22; 
nShot        = 419; 
% The FoV can be set as 480 here (low-resolution), 
% since we want to estimate coil sensitivity
Matrix_size = 480;
acquisitonFoV          = [240, 240, 240];
% This number (nShotOff) has to be adapted based on the observation of the 
% steady-state graph
nShotOff     = 14; 
N_u          = [48, 48, 48]; % the size of k-space
reconFov = 240;
dK_u         = [1, 1, 1]./reconFov; % interval between k-space

if subject_num == 1
    nCh_array    = 42;
elseif subject_num == 2
    nCh_array    = 44;
elseif subject_num == 3
    nCh_array = 34;
else
    nCh_array = 42;
end

nCh_body     = 2; 
% We will need to ask for a predefined format. Hence we need to read the data outside and
% pass the y_body, t has to be computed by the user.  
% We use trajectory in using the phisical dimentions without convention [-0.5,0.5]
% You have to see how much data needs to be discarded by looking at the
% graph, you don't want to keep data if you are not in steady-state.
%%
[y_body, t, ve] = bmCoilSense_nonCart_dataFromTwix( bodyCoilFile, ...
                                                    N_u, ...
                                                    N, ...
                                                    nSeg, ...
                                                    nShot, ...
                                                    nCh_body, ...
                                                    reconFov, ...
                                                    nShotOff);
% same for arraycoil 
y_array         = bmCoilSense_nonCart_dataFromTwix( arrayCoilFile, ...
                                                    N_u, ...
                                                    N, ...
                                                    nSeg, ...
                                                    nShot, ...
                                                    nCh_array, ...
                                                    reconFov, ...
                                                    nShotOff);
% compute the gridding matrices (Gn = approximation of inverse, Gu = Forward,
% Gut = transposed of Gu) Gn and Gut are both backward
[Gn, Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u); 
disp('Gn, Gu, Gut are computed')
%
% You need to Reassign the xmin, xmax & ymin, ymax & zmin, zmax 
% To do it you need to run the function below (bmCoilSense_nonCart_mask)
% control + E: to change the tresholds
% shift + E: to set the constast chosen

% Box excluding coordinates
x_min = 1; 
x_max = 40;

y_min = 9; 
y_max = 41;

z_min = 4; 
z_max = 48;

% Two thresholds
th_RMS = 14; 
th_MIP = 10; 

close_size = []; 
open_size  = []; 
% Mask computation, the two thresholds to exclude artifacts from the region
% where there is no signal, like air in the lungs: 
% 1 for the Root mean square
% 1 for the Maximum intensity projection

% Make some work to try to automate it. Define two scripts.
% (Automatic and Advanced)
m = bmCoilSense_nonCart_mask(   y_body, Gn, ...
                                x_min, x_max, ...
                                y_min, y_max, ...
                                z_min, z_max, ...
                                th_RMS, th_MIP, ...
                                close_size, ...
                                open_size, ...
                                false);
% Reset the parameters in the cell above and rerun this cell.
% Figure 1:
    % X axis: the normalized intensity 
    % Y axis: the percentage of pixels whose internsities are above the normalized intensity
    % The thresholds should be at the elbow of the curve (most of the pixels)

%%

% Select one body coil and compute its sensitivity
[y_ref, C_ref] = bmCoilSense_nonCart_ref(y_body, Gn, m, []); 



% Estimate the coil sensitivity of each surface coil using one body coil
% image as reference image C_c = (X_c./x_ref)
C_array_prime = bmCoilSense_nonCart_primary(y_array, y_ref, C_ref, Gn, ve, m);


% Do a recon, predending the selected body coil is one channel among the
% others, and optimize the coil sensitivity estimate by alternating steps
% Of gradient descent (X,C)
nIter = 5; 
[C, x] = bmCoilSense_nonCart_secondary(y_array, C_array_prime, y_ref, C_ref, Gn, Gu, Gut, ve, nIter, false); 
close all;
%% Save C into the folder
if subject_num == 1
    saveCDirList = {'/Sub001/T1_LIBRE_Binning/C/','/Sub001/T1_LIBRE_woBinning/C/'};
elseif subject_num == 2
    saveCDirList = {'/Sub002/T1_LIBRE_Binning/C/','/Sub002/T1_LIBRE_woBinning/C/'};
elseif subject_num == 3
    saveCDirList = {'/Sub003/T1_LIBRE_Binning/C/','/Sub003/T1_LIBRE_woBinning/C/'};
else
    saveCDirList = {'/Sub004/T1_LIBRE_Binning/C/','/Sub004/T1_LIBRE_woBinning/C/'};
end

for idx = 1:2
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
%%
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
%% Generate mask


    % TimeStamp difference, see if it matches the ET data
% T1-weighted
% name psychopy and ET data: JB1
% 
% duration of the protocol: 720sec (12 min)
% 
% Our protocol (p178 362s + 15.03s + 357s)
% 
% differences between the two sequence: 15.03 sec
% 
% ET calibration: 5 points
% 
% LIBRE Voxel size: 0.5mm * 0.5mm * 0.5mm
% 
% VIBE Voxel size: 0.4mm * 0.4mm * 0.4mm


% In protocol, we set the T1 LIBRE duration: 362(sec)*1000(freq)
% Sub001: However, here the duration of the rawdata is: 374.5575sec with data points:45210
% Sub002: However, here the duration of the rawdata is: 374380 ms with data points:45210
% Sub003: 334237.5ms with data points:40340
% Sub004: However, here the duration of the rawdata is: 374565 ms with data points:45210
%% Some information given by the Eyetracker data
disp('If you have not got mask from ET prepared, go to ET analysis.')
% No need to select if using standard seq.
% x direction and y direction
% nLine 
    
%--------------------------------------------------------------------------
% Assign a bin to each line (X)
%--------------------------------------------------------------------------
th = 7; 
nShotOff = 14; 
nSeg = 22; 
cMask = eyeGenerateBinning(datasetDir,nShotOff, nSeg, th, ETDir, true);

%% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    

%save(fullfile(param.savedir,'cMask.mat'),'param');
%disp(['param is saved here:', param.savedir, '\cMask.mat'])
if subject_num == 1
    otherDir = [reconDir, '/Sub001/T1_LIBRE_Binning/other/'];
elseif subject_num == 2
    otherDir = [reconDir, '/Sub002/T1_LIBRE_Binning/other/'];
elseif subject_num == 3
    otherDir = [reconDir, '/Sub003/T1_LIBRE_Binning/other/'];
else
    otherDir = [reconDir, '/Sub004/T1_LIBRE_Binning/other/'];
end
% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

% CMaskFilePath = [otherDir, sprintf('cMask_th%d.mat', th)];
CMaskFilePath = [otherDir, 'cMask_4fr.mat'];
% Save the CMask to the .mat file
save(CMaskFilePath, 'cMask');
disp('cMask has been saved here:')
disp(CMaskFilePath)

%%
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% Mitosius
% clearvars -except reconDir datasetDir nbins nShotOff th CMaskFilePath otherDir subject_num reconFov

% path to the coil sensitivity C.mat
if subject_num == 1
    CfilePath = fullfile(reconDir, '/Sub001/T1_LIBRE_Binning/C/C.mat');
elseif subject_num == 2
    CfilePath = fullfile(reconDir, '/Sub002/T1_LIBRE_Binning/C/C.mat');
elseif subject_num == 3
    CfilePath = fullfile(reconDir, '/Sub003/T1_LIBRE_Binning/C/C.mat');
else
    CfilePath = fullfile(reconDir, '/Sub004/T1_LIBRE_Binning/C/C.mat');
end
% th=7
%subject_num = 1
if subject_num == 1
    CMaskFilePath = [reconDir, sprintf('/Sub001/T1_LIBRE_Binning/other/cMask_th%d.mat', th)];
% path to the rawdatafile (in this case Siemens raw data)
    data_filename = '/meas_MID00453_FID57919_BEAT_LIBREon_eye.dat';
elseif subject_num == 2
    CMaskFilePath = [reconDir, sprintf('/Sub002/T1_LIBRE_Binning/other/cMask_th%d.mat', th)];
% path to the rawdatafile (in this case Siemens raw data)
    data_filename = '/meas_MID00342_FID56821_BEAT_LIBREon_eye.dat';
elseif subject_num == 3
    CMaskFilePath = [reconDir, sprintf('/Sub003/T1_LIBRE_Binning/other/cMask_th%d.mat', th)];
% path to the rawdatafile (in this case Siemens raw data)
    data_filename = '/meas_MID00297_FID57763_BEAT_LIBREon_eye.dat';
else
    CMaskFilePath = [reconDir, sprintf('/Sub004/T1_LIBRE_Binning/other/cMask_th%d.mat', th)];
% path to the rawdatafile (in this case Siemens raw data)
    data_filename = '/meas_MID00333_FID57799_BEAT_LIBREon_eye.dat';
end


%
% Initialize and fill in the parameters: This in theory can be automated;
% However Bastien told us that it can lead to errors. In addition we should
% be indipendent from the specific file format used.
p = bmMriAcquisitionParam([]); 
p.name            = [];
p.mainFile_name   = data_filename;

p.imDim           = 3;
p.N     = 480;  
p.nSeg  = 22;  
p.nShot = 2055;  
p.nLine = 45210;  
p.nPar  = 1;  

p.nLine           = double([]);
p.nPt             = double([]);
p.raw_N_u         = [Matrix_size, Matrix_size, Matrix_size];
p.raw_dK_u        = [1, 1, 1]./reconFov;

if subject_num == 1
    p.nCh   = 42;%note
elseif subject_num == 2
    p.nCh   = 44;%note  
elseif subject_num == 3
    p.nCh   = 34;
else
    p.nCh   = 42;
end
p.nEcho = 1; 

p.selfNav_flag    = true;
% This was estimated in the coil sensitivity computation， the first 10
% shots are eliminated.
% nShotOff = 10;
p.nShot_off       = nShotOff; 
p.roosk_flag      = false;
% This is the full FOV not the half FOV
p.FoV             = [reconFov, reconFov, reconFov];
% This sets the trajectory used
p.traj_type       = 'full_radial3_phylotaxis';

% Fill in missing parameters that can be deduced from existing ones.
p.refresh; 

%% read rawdata
f = [datasetDir, data_filename]; 
% % Display infos
bmTwix_info(f); 
% read raw data
myTwix = bmTwix(f);
y_tot   = bmTwix_data(myTwix, p);
% compute trajectory points. This function is really wird. ASK BASTIEN.
t_tot   = bmTraj(p); 
% compute volume elements
ve_tot  = bmVolumeElement(t_tot, 'voronoi_full_radial3' ); 


N_u     = [Matrix_size, Matrix_size, Matrix_size];
n_u     = [Matrix_size, Matrix_size, Matrix_size];
dK_u    = [1, 1, 1]./reconFov;

% Load the coil sensitivity previously measured
load(CfilePath);

C = bmImResize(C, [48, 48, 48], N_u); 

%% Normalization (probably to converge better)
% Note you can normalize the rawdata and the image will be normalized
% This is because the Fourier transform is linear
% F(f(.)/a) =  F(f(.))/a
if Matrix_size >240
    normalization = false;
end
if normalization
    x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, n_u, dK_u); 
    %
    bmImage(x_tot)
    temp_im = getimage(gca); 
    
    bmImage(temp_im); 
    temp_roi = roipoly; 
    normalize_val = mean(temp_im(temp_roi(:))); 
    % The normalize_val is super small, it is 5e-10, very small
    % The value of one complex point is like: -0.0396 - 0.1162i
    y_tot(1,1,123)
end
% only once !!!!
if real(y_tot)<1
    if normalization
        y_tot = y_tot/normalize_val; 
        y_tot(1,1,123)
    else
        y_tot = y_tot/(5e-9); 
        y_tot(1,1,123)
    end
end
%% Before running this cell, make sure the cMask is well-prepared.
% Load the masked coil sensitivity 

load(CMaskFilePath); 
disp(CMaskFilePath);
disp('is loaded!')

% Eleminate the first segment of all the spokes for accuracies
cMask = reshape(cMask, [nbins, 22, 2055]); 
cMask(:, 1, :) = []; 

cMask(:, :, 1:p.nShot_off) = []; 
cMask = bmPointReshape(cMask); 

%%
% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot, t_tot, cMask); 
y = bmPermuteToCol(y); 

ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

%% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction
if subject_num == 1
    mDir = [reconDir, '/Sub001/T1_LIBRE_Binning/mitosius/th7'];
elseif subject_num == 2
    mDir = [reconDir, '/Sub002/T1_LIBRE_Binning/mitosius/th7'];
elseif subject_num == 3
    mDir = [reconDir, '/Sub003/T1_LIBRE_Binning/mitosius/th7'];
else
    mDir = [reconDir, '/Sub004/T1_LIBRE_Binning/mitosius/th7'];
end
bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!');
disp(mDir);
%% Now generate the cMask for wo Binning 
% And rerun mitusius
%--------------------------------
disp('Generating mask for wo Binning')
binMaskMatrix_woBinning = ones([p.nLine ,nbins]);
nMeasuresOff = nShotOff * p.nSeg;
for k=1:p.nLine
    if k<= nMeasuresOff
        binMaskMatrix_woBinning(k, 1) = 0;
    end
    if mod(k,22) == 1
        binMaskMatrix_woBinning(k, 1) = 0;
    end
end
sum_binning = sum(binMaskMatrix_woBinning);
disp(['wo Binning, preserved #line',num2str(sum_binning)])
cMask_woBinning = logical(binMaskMatrix_woBinning)';
if subject_num == 1
    otherDir = [reconDir, '/Sub001/T1_LIBRE_woBinning/other/'];
elseif subject_num == 2
    otherDir = [reconDir, '/Sub002/T1_LIBRE_woBinning/other/'];
elseif subject_num == 3
    otherDir = [reconDir, '/Sub003/T1_LIBRE_woBinning/other/'];
else
    otherDir = [reconDir, '/Sub004/T1_LIBRE_woBinning/other/'];
end
% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end
CMaskwoBinningFilePath = [otherDir, 'cMask_woBinning.mat'];
save(CMaskwoBinningFilePath, 'cMask_woBinning');
disp('cMask_woBinning has been saved here:')
disp(CMaskwoBinningFilePath)

% Eleminate the first segment of all the spokes for accuracies
cMask_woBinning = reshape(cMask_woBinning, [nbins, 22, 2055]); 
cMask_woBinning(:, 1, :) = []; 

cMask_woBinning(:, :, 1:p.nShot_off) = []; 
cMask_woBinning = bmPointReshape(cMask_woBinning); 

% Run the mitosis function and compute volume elements
[y, t] = bmMitosis(y_tot, t_tot, cMask_woBinning); 
y = bmPermuteToCol(y); 

ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 
if subject_num == 1
    mDir = [reconDir, '/Sub001/T1_LIBRE_woBinning/mitosius'];
elseif subject_num == 2
    mDir = [reconDir, '/Sub002/T1_LIBRE_woBinning/mitosius'];
elseif subject_num == 3
    mDir = [reconDir, '/Sub003/T1_LIBRE_woBinning/mitosius'];
else
    mDir = [reconDir, '/Sub004/T1_LIBRE_woBinning/mitosius'];
end
bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)

%%
% scp -r "/media/sinf/1,0 TB Disk/240922_480/Sub001/T1_LIBRE_Binning/mitosius" yi9226@hpc1-login.chuv.ch:/data/line/MREye_yiwei/240922_480/Sub001/T1_LIBRE_Binning
% scp -r "/media/sinf/1,0 TB Disk/240922_480/Sub001/T1_LIBRE_woBinning/mitosius" yi9226@hpc1-login.chuv.ch:/data/line/MREye_yiwei/240922_480/Sub001/T1_LIBRE_woBinning
% scp -r "/media/sinf/1,0 TB Disk/240922_480/Sub001/T1_LIBRE_Binning/C" yi9226@hpc1-login.chuv.ch:/data/line/MREye_yiwei/240922_480/Sub001/T1_LIBRE_Binning
% scp -r "/media/sinf/1,0 TB Disk/240922_480/Sub001/T1_LIBRE_woBinning/C" yi9226@hpc1-login.chuv.ch:/data/line/MREye_yiwei/240922_480/Sub001/T1_LIBRE_Binning
disp(sum(cMask_woBinning))