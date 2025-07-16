clc; clear all;
debug = false;

Matrix_size = 240;
reconFov = 240;
region_idx = 0;
% 0:up 1:down 2:left 3:right 4:center mask

% change woBinning/Binning, change th75->th0, change eyeMask name!
%%
subject_num=1;
reconDir = '/home/debi/jaime/tmp/250613_JB/';

if subject_num == 1
    datasetDir = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-01/rawdata/';
    ETDir      = '/home/debi/jaime/repos/MR-EyeTrack/data/pilot/sub-01/masks_1206/';
end


if subject_num == 1
    bodyCoilFile     = [datasetDir, '/meas_MID2400614_FID182868_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00615_FID182869_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00605_FID182859_BEAT_LIBREon_eye_(23_09_24).dat'];
    
elseif subject_num == 2
    bodyCoilFile     = [datasetDir, '/meas_MID00589_FID182843_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00590_FID182844_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00580_FID182834_BEAT_LIBREon_eye_(23_09_24).dat'];
elseif subject_num == 3
    bodyCoilFile = [datasetDir, '/meas_MID00563_FID182817_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile = [datasetDir, '/meas_MID00564_FID182818_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00554_FID182808_BEAT_LIBREon_eye_(23_09_24).dat'];

end


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


if subject_num == 1
    saveCDirList = {'/Sub001/T1_LIBRE_Binning/C/','/Sub001/T1_LIBRE_woBinning/C/'};
elseif subject_num == 2
    saveCDirList = {'/Sub002/T1_LIBRE_Binning/C/','/Sub002/T1_LIBRE_woBinning/C/'};
elseif subject_num == 3
    saveCDirList = {'/Sub003/T1_LIBRE_Binning/C/','/Sub003/T1_LIBRE_woBinning/C/'};
else
    saveCDirList = {'/Sub004/T1_LIBRE_Binning/C/','/Sub004/T1_LIBRE_woBinning/C/'};
end



%
myTwix = bmTwix(measureFile); 

%
p = bmMriAcquisitionParam([]); 
p.name            = [];
p.mainFile_name   = measureFile;

p.imDim           = 3;
p.N     = 480;  
p.nSeg  = 22;  
p.nShot = 3723;  
p.nLine = 81906;   
p.nPar  = 1;  

p.nLine           = double([]);
p.nPt             = double([]);
p.raw_N_u         = [Matrix_size, Matrix_size, Matrix_size];
p.raw_dK_u        = [1, 1, 1]./reconFov;


 
if subject_num == 1
    p.nCh    = 44;
elseif subject_num == 2
    p.nCh    = 52;
elseif subject_num == 3
    p.nCh = 52;
else
    p.nCh = 999;
end

p.nEcho = 1; 

p.selfNav_flag    = true;
% This was estimated in the coil sensitivity computation， the first 10
% shots are eliminated.
p.nShot_off       = 15; 
p.roosk_flag      = false;
% This is the full FOV not the half FOV
p.FoV             = [reconFov, reconFov, reconFov];
% This sets the trajectory used
p.traj_type       = 'full_radial3_phylotaxis';

% Fill in missing parameters that can be deduced from existing ones.
p.refresh; 

% read rawdata
y_tot   = bmTwix_data(myTwix, p);

% compute trajectory points. This function is really wird. ASK BASTIEN.
t_tot   = bmTraj(p); 
% compute volume elements
ve_tot  = bmVolumeElement(t_tot, 'voronoi_full_radial3' ); 


N_u     = [Matrix_size, Matrix_size, Matrix_size];
n_u     = [Matrix_size, Matrix_size, Matrix_size];
dK_u    = [1, 1, 1]./reconFov;

% Load the coil sensitivity previously measured
saveCDir     = [reconDir,saveCDirList{1}];
CfileName = 'C.mat';
CfilePath = fullfile(saveCDir, CfileName);

load(CfilePath);

C = bmImResize(C, [48, 48, 48], N_u); 

% Normalization (probably to convege better)
if Matrix_size >240
    normalization = false;
else 
    normalization = true;
end
%
if normalization
    x_tot = bmMathilda(y_tot, t_tot, ve_tot, C, N_u, n_u, dK_u); 
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
        y_tot = y_tot/(5e-9); 
        y_tot(1,1,123)
    end
end
close all;

%% Before running this cell, make sure the Mask is well-prepared.
% Load the masked coil sensitivity 

% if woBinning
% MaskFilePath = [otherDir, 'eMask_woBin.mat'];
% if withBinning

for region_idx = 0:3

    if subject_num == 1
        mDir = [reconDir, '/Sub001/T1_LIBRE_Binning/mitosius/mask_', num2str(region_idx)];
    elseif subject_num == 2
        mDir = [reconDir, '/Sub002/T1_LIBRE_Binning/mitosius/mask_', num2str(region_idx)];
    elseif subject_num == 3
        mDir = [reconDir, '/Sub003/T1_LIBRE_Binning/mitosius/mask_', num2str(region_idx)];
    else
        mDir = [reconDir, '/Sub004/T1_LIBRE_Binning/mitosius/mask_', num2str(region_idx)];
    end

    th_ratio = 3/4;
    eMaskFilePath = [otherDir,sprintf('eMask_th%.2f_region%d.mat', th_ratio, region_idx)];
    %
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
    
    
    % Run the mitosis function and compute volume elements
    
    [y, t] = bmMitosis(y_tot, t_tot, eyeMask); 
    y = bmPermuteToCol(y); 
    ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

    % Save all the resulting datastructures on the disk. You are now ready
    % to run your reconstruction
    
    bmMitosius_create(mDir, y, t, ve); 
    disp('Mitosius files are saved!')
    disp(mDir)

end




















