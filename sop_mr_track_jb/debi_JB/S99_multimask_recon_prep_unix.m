clear all;
clc;
%%
% Coil sensitivity -> binning mask cMask -> Mitosius
% ------------------------------------------------------------------
% ------------------------------------------------------------------

% ------------------------------------------------------------------
% ------------------------------------------------------------------
%% Initialize the directories and acquire the Coil
subject_num = 1;
%%
reconDir = '/media/sinf/1,0 TB Disk/240912_recon/';
if subject_num == 1
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/RawData';
elseif subject_num == 2
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/MR_EYE_Subj02/RawData';
elseif subject_num == 3
    datasetDir = 'C:\yiwei\1_Pilot_MREye_Data\Sub001\230928_anatomical_MREYE_study\MR_EYE_Subj03\RawData';
else
    datasetDir = 'C:\yiwei\1_Pilot_MREye_Data\Sub001\230928_anatomical_MREYE_study\MR_EYE_Subj04\RawData';
end

if subject_num == 1
    bodyCoilFile     = [datasetDir, '/meas_MID00469_FID57935_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00470_FID57936_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00453_FID57919_BEAT_LIBREon_eye.dat'];
else
    bodyCoilFile     = [datasetDir, '/meas_MID00357_FID56836_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00358_FID56837_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00342_FID56821_BEAT_LIBREon_eye.dat'];
end
%% Start coil estimation
disp('-------------------------------------')
disp('Start Coil Estimation')
disp('-------------------------------------')
% Read metadata from the twix
bmTwix_info(bodyCoilFile)
bmTwix_info(arrayCoilFile)
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
acquisitonFoV          = [240, 240, 240];
% This number (nShotOff) has to be adapted based on the observation of the 
% steady-state graph
nShotOff     = 14; 
N_u          = [48, 48, 48]; % the size of k-space
reconFov = 240;
dK_u         = [1, 1, 1]./reconFov; % interval between k-space

if subject_num == 1
    nCh_array    = 42;
else
    nCh_array    = 44;%note subject 2 here nch=44 
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
%%
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

%% Make some work to try to automate it. Define two scripts.
% (Automatic and Advanced)
m = bmCoilSense_nonCart_mask(   y_body, Gn, ...
                                x_min, x_max, ...
                                y_min, y_max, ...
                                z_min, z_max, ...
                                th_RMS, th_MIP, ...
                                close_size, ...
                                open_size, ...
                                true);
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
[C, x] = bmCoilSense_nonCart_secondary(y_array, C_array_prime, y_ref, C_ref, Gn, Gu, Gut, ve, nIter, true); 

%% Save C into the folder
if subject_num == 1
    saveCDirList = {'/Sub001/T1_LIBRE_Binning/C/','/Sub001/T1_LIBRE_woBinning/C/'};
else
    saveCDirList = {'/Sub002/T1_LIBRE_Binning/C/','/Sub002/T1_LIBRE_woBinning/C/'};
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
%--------------------------------------------------------------------------
% custom change to import mapVBVD: not sure its correct, no time to test it.
% Just need to be sure you are able to acces the mapVBVD_JH function for example.
%--------------------------------------------------------------------------
addpath('/home/sinf/yiwei/Recon_fork/MR_EYE_RawData_Binning_fork/code/MATLAB/ReadRawDataSiemens/mapVBVD')
param.basedir = datasetDir;
basedir = param.basedir;
param.batchParam = [];

%%
% Prompt user to select Siemens Raw Data file
%--------------------------------------------------------------------------

    [rawDataName, rawDataDir, ~] = uigetfile( ...
    { '*.dat','Siemens raw data file (*.dat)'}, ...
       'Pick a file', ...
       'MultiSelect', 'off',basedir);

    if rawDataName == 0
        warning('No file selected');
        return;
    end

    filepathRawData = fullfile(rawDataDir, rawDataName);
%%
%--------------------------------------------------------------------------
% Initialize saving directory
%--------------------------------------------------------------------------    
  % Create saving directory
    [~,name,~] = fileparts(rawDataName);

  % Add number of bin to saving directory name
    name = sprintf('BinningEYE_%s',name);
    if subject_num == 1
    saveBinningDir = [reconDir, '/Sub001/T1_LIBRE_Binning/'];
    elseif subject_num == 2
    saveBinningDir = [reconDir, '/Sub002/T1_LIBRE_Binning/'];
    end
    param.savedir = fullfile(saveBinningDir,name);
    disp(param.savedir) 
  % Create directory
    if exist(param.savedir,'dir') ~= 7
        mkdir(param.savedir);
    end
    disp(param.savedir)
    if ~exist(param.savedir, 'dir')
        mkdir(param.savedir);
    end
%%
%--------------------------------------------------------------------------
% Prompt user to enter the number of desired bins
%-------------------------------------------------------------------------- 
% How to select the number of bins?
    prompt        = {'Enter the number of bins'};
    name          = '#Bins';
    numlines      = 1;
    defaultanswer = {'1'};
    answer        = inputdlg(prompt,name,numlines,...
                           defaultanswer);
    
  % Check if the user selected cancel
    if isempty(answer)
        warning('The user selected cancel');
        return;
    end
    
  % Convert string answer to number
    nbins = str2double(answer{1});
    param.nBins = nbins;

   
%--------------------------------------------------------------------------
% Read the PMUTimeStamp triggered from external trigger
%--------------------------------------------------------------------------  
    costTime = 2.5;
    param.batchParam.rawDataName    = rawDataName;
    param.batchParam.rawDataDir     = rawDataDir;

    [ twix_obj, param ] = dataSelectionAndLoading( basedir, param );
    % Load the PMUTime and TimeStamp, shift the TimeStamp to the beginning
    % of 0
    PMUTimeStamp    = double( twix_obj.image.pmutime );
    TimeStamp       = double( twix_obj.image.timestamp );

    %
    TimeStamp       = TimeStamp - min(TimeStamp);
    % Do not forget to scale the times by costTime (Setting from Siemens)
    PMUTimeStamp_ms = PMUTimeStamp * costTime;
    PMUTimeStamp_s  = PMUTimeStamp_ms / 1000;
    TimeStamp_ms    = TimeStamp * costTime;
    TimeStamp_s     = TimeStamp_ms / 1000;
    % Save all the param to struct "param"
    param.PMUTimeStamp_ms   = PMUTimeStamp_ms;
    param.PMUTimeStamp_s    = PMUTimeStamp_s;
    param.TimeStamp_ms      = TimeStamp_ms;
    param.TimeStamp_s       = TimeStamp_s;

%% Inspect the visualization of PMUTime and TimeStamp
inspect = false;
if inspect
    figure;
    plot(param.TimeStamp_ms, param.PMUTimeStamp_ms);
    % xlim([0 1e4]);  % Set x-axis limits
    % ylim([0 180]); % Set y-axis limits
    xlabel('TimeStamp (ms)');
    ylabel('PMU TimeStamp (ms)');
    title('Time vs. PMU Stamp');
end
    %% TimeStamp difference, see if it matches the ET data
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
Timediff = TimeStamp_ms(end) - TimeStamp_ms(1);

disp('In protocol, we set the T1 LIBRE duration: 362(sec)*1000(freq)')
disp('And we set the T1 VIBE duration: 357(sec)*1000(freq)')
disp(['However, here the duration of the rawdata is: ', num2str(Timediff), ...
    ' ms with data points:', num2str(length(TimeStamp_s)) ]);
% In protocol, we set the T1 LIBRE duration: 362(sec)*1000(freq)
% Sub001: However, here the duration of the rawdata is: 374.5575sec with data points:45210
% Sub002: However, here the duration of the rawdata is: 374380 ms with data points:45210
%% Some information given by the Eyetracker data
disp('If you have not got mask from ET prepared, go to ET analysis.')
% No need to select if using standard seq.
% x direction and y direction
% nLine 
%%
% Prompt user to select generated mask file from ET data (.mat)
%--------------------------------------------------------------------------

[maskDataName, maskDataDir, ~] = uigetfile( ...
{ '*.mat','Generated mask data file (*.mat)'}, ...
   'Pick a file', ...
   'MultiSelect', 'off', reconDir);

if maskDataName == 0
    warning('No mask file selected');
    return;
end

filepathMaskData = fullfile(maskDataDir, maskDataName);
mask_method_1 = load(filepathMaskData);
if isstruct(mask_method_1)
    mask_method_1 = mask_method_1.array;
end

%%     
%--------------------------------------------------------------------------
% Assign a bin to each line (X)
%--------------------------------------------------------------------------
    NLin    = length(param.PMUTimeStamp_ms);
    binMask = cell(nbins,1);
    
    % Here nbins is 1
    binMaskMatrix = zeros([NLin,nbins]);
    % number of timestamps, number of bins

    % Set the window width and threshold
    WinWidth = round(length(mask_method_1)/NLin);
    th = WinWidth-1;

    %% Initialize
    timeSegMinus1 = -1;
    timeSeg = 0;
    timeSegPlus1 = 0;
    % nShotOff=10
    % nSeg=22
    nMeasuresOff = nShotOff*nSeg;
    for k = 1:NLin
        % Eliminate the first nShotOff
        if k<= nMeasuresOff
             binMaskMatrix( k, 1 ) = 0;
        else
            if k == NLin
                timeSegMinus1 = TimeStamp_ms(k-1);
                timeSeg = TimeStamp_ms(k);
                timeSegPlus1 = TimeStamp_ms(k);
            else
                timeSegMinus1 = TimeStamp_ms(k-1);
                timeSeg = TimeStamp_ms(k);
                timeSegPlus1 = TimeStamp_ms(k+1);
            end
            
    
            win_lower = ceil(1/2*(timeSegMinus1+timeSeg));
            if win_lower < 1
                win_lower = 1;
            end
            win_upper = floor(1/2*(timeSegPlus1+timeSeg));
            if win_upper > numel(mask_method_1)
                win_upper = numel(mask_method_1);
            end
    
            window_data = mask_method_1(win_lower:win_upper);
            if mod(k,100)==1
                disp(window_data);
            end
            % Count the number of True values
            true_count = sum(window_data);
            if true_count >= th
                binMaskMatrix( k, 1 ) = 1;
            end
        end

        if mod(k, 22)  == 1
            binMaskMatrix( k, 1 ) = 0;
        end
    end
    
    
        
    for k = 1:nbins
        binMask{k} = binMaskMatrix(:,k);
    end
    sum_binning = sum(binMaskMatrix);
    disp(['threshold: ', num2str(th)]);
    disp(['with Binning, preserved #line: ',num2str(sum_binning)])
    param.binMask = binMask;
    % th: 6 preserved:39303/45210
    % th: 7 preserved: 39135/45210
    % noBinning: 42945/45210

%% Saving data and Convert to Monalisa format
%--------------------------------------------------------------------------    

%save(fullfile(param.savedir,'cMask.mat'),'param');
%disp(['param is saved here:', param.savedir, '\cMask.mat'])
if subject_num == 1
    otherDir = [reconDir, '/Sub001/T1_LIBRE_Binning/other/'];
else
    otherDir = [reconDir, '/Sub002/T1_LIBRE_Binning/other/'];
end
% Check if the directory exists
if ~isfolder(otherDir)
    % If it doesn't exist, create it
    mkdir(otherDir);
    disp(['Directory created: ', otherDir]);
else
    disp(['Directory already exists: ', otherDir]);
end

cMask = logical(binMaskMatrix)';
CMaskFilePath = [otherDir, sprintf('cMask_th%d.mat', th)];
% Save the CMask to the .mat file
save(CMaskFilePath, 'cMask');
disp('cMask has been saved here:')
disp(CMaskFilePath)

%%
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
% ---------------------------------------------------------------------------------
clc;
clearvars -except reconDir datasetDir nbins nShotOff th CMaskFilePath otherDir subject_num reconFov

% path to the coil sensitivity C.mat
if subject_num == 1
    CfilePath = fullfile(reconDir, '/Sub001/T1_LIBRE_Binning/C/C.mat');
else
    CfilePath = fullfile(reconDir, '/Sub002/T1_LIBRE_Binning/C/C.mat');
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
end

f = [datasetDir, data_filename]; 
% % Display infos
bmTwix_info(f); 
% read raw data
myTwix = bmTwix(f);
%%
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
Matrix_size = 240;
p.nLine           = double([]);
p.nPt             = double([]);
p.raw_N_u         = [Matrix_size, Matrix_size, Matrix_size];
p.raw_dK_u        = [1, 1, 1]./reconFov;

if subject_num == 1
    p.nCh   = 42;%note
else
    p.nCh   = 44;%note  
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
x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, n_u, dK_u); 
%%
bmImage(x_tot)
temp_im = getimage(gca); 

bmImage(temp_im); 
temp_roi = roipoly; 
normalize_val = mean(temp_im(temp_roi(:))); 
% The normalize_val is super small, it is 5e-10, very small
% The value of one complex point is like: -0.0396 - 0.1162i
y_tot(1,1,123)
%% only once !!!!
if real(y_tot)<1
    y_tot = y_tot/normalize_val; 
    y_tot(1,1,123)
else
    y_tot(1,1,123)
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


% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot, t_tot, cMask); 
y = bmPermuteToCol(y); 

ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

%% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction
if subject_num == 1
    mDir = [reconDir, '/Sub001/T1_LIBRE_Binning/mitosius/th7'];
else
    mDir = [reconDir, '/Sub002/T1_LIBRE_Binning/mitosius/th7'];
end
bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!');
disp(mDir);
%% Now generate the cMask for wo Binning 
% And rerun mitusius
%--------------------------------
disp('Generating mask for wo Binning')
binMaskMatrix = ones([p.nLine ,nbins]);
nMeasuresOff = nShotOff * p.nSeg;
for k=1:p.nLine
    if k<= nMeasuresOff
        binMaskMatrix(k, 1) = 0;
    end
    if mod(k,22) == 1
        binMaskMatrix(k, 1) = 0;
    end
end
sum_binning = sum(binMaskMatrix);
disp(['wo Binning, preserved #line',num2str(sum_binning)])
cMask_woBinning = logical(binMaskMatrix)';
if subject_num == 1
    otherDir = [reconDir, '/Sub001/T1_LIBRE_woBinning/other/'];
else
    otherDir = [reconDir, '/Sub002/T1_LIBRE_woBinning/other/'];
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
else
    mDir = [reconDir, '/Sub002/T1_LIBRE_woBinning/mitosius'];
end
bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)

%%

disp(sum(cMask_woBinning))