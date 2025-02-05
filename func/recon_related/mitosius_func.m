function mitosius_func(subject_num, c_path, cMask_path, mDir)
debug = false;
% subject_num = 4;
% reconDir = '/media/sinf/1,0 TB Disk/240922_480/';
if subject_num == 1
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub001/230928_anatomical_MREYE_study/MR_EYE_Subj01/RawData';
elseif subject_num == 2
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub002/230926_anatomical_MREYE_study/MR_EYE_Subj02/RawData';
elseif subject_num == 3
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub003/230928_anatomical_MREYE_study/MR_EYE_Subj03/RawData';
else
    datasetDir = '/media/sinf/1,0 TB Disk/1_Pilot_MREye_Data/Sub004/230923_anatomical_MREYE_study/MR_EYE_Subj04/RawData';
end

if subject_num == 1
    bodyCoilFile     = [datasetDir, '/meas_MID00469_FID57935_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00470_FID57936_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00453_FID57919_BEAT_LIBREon_eye.dat'];
    rawDataName = 'meas_MID00453_FID57919_BEAT_LIBREon_eye.dat';
elseif subject_num == 2
    bodyCoilFile     = [datasetDir, '/meas_MID00357_FID56836_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile    = [datasetDir, '/meas_MID00358_FID56837_BEAT_LIBREon_eye_HC_BC.dat'];

    measureFile = [datasetDir, '/meas_MID00342_FID56821_BEAT_LIBREon_eye.dat'];
    rawDataName = 'meas_MID00342_FID56821_BEAT_LIBREon_eye.dat';
elseif subject_num == 3
    bodyCoilFile = [datasetDir, '/meas_MID00312_FID57778_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile = [datasetDir, '/meas_MID00313_FID57779_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00297_FID57763_BEAT_LIBREon_eye.dat'];
    rawDataName = 'meas_MID00297_FID57763_BEAT_LIBREon_eye.dat';
else
    bodyCoilFile = [datasetDir, '/meas_MID00349_FID57815_BEAT_LIBREon_eye_BC_BC.dat'];
    arrayCoilFile = [datasetDir, '/meas_MID00350_FID57816_BEAT_LIBREon_eye_HC_BC.dat'];
    measureFile = [datasetDir, '/meas_MID00333_FID57799_BEAT_LIBREon_eye.dat'];
    rawDataName = 'meas_MID00333_FID57799_BEAT_LIBREon_eye.dat';
end

Matrix_size = 480;

reconFov = 240;
% c_path = [reconDir, '/Sub004/T1_LIBRE_Binning/C/C.mat'];
% cMask_path = [reconDir, '/Sub004/T1_LIBRE_woBinning/other/cMask_woBinning_sampled_12912.mat'];
% mDir = '/media/sinf/1,0 TB Disk/240922_480/Sub004/T1_LIBRE_woBinning/mitosius/woBinning_sampled_12912/';
if ~exist(mDir, 'dir')
    mkdir(mDir);
end
    
%
% Display infos
% bmTwix_info(measureFile); 
% read raw data
myTwix = bmTwix(measureFile); 
%
% Initialize and fill in the parameters: This in theory can be automated;
% However Bastien told us that it can lead to errors. In addition we should
% be indipendent from the specific file format used.
p = bmMriAcquisitionParam([]); 
p.name            = [];
p.mainFile_name   = measureFile;

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
    p.nCh = 34;
else
    p.nCh = 42;
end

p.nEcho = 1; 

p.selfNav_flag    = true;
% This was estimated in the coil sensitivity computation， the first 10
% shots are eliminated.
p.nShot_off       = 14; 
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
load(c_path);

C = bmImResize(C, [48, 48, 48], N_u); 

% Normalization (probably to convege better)
if Matrix_size >240
    normalization = false;
else 
    normalization = true;
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


% Before running this cell, make sure the cMask is well-prepared.
% Load the masked coil sensitivity 
cMask = load(cMask_path); 
fields = fieldnames(cMask);  % Get the field names
firstField = fields{1};  % Get the first field name
cMask = cMask.(firstField);  % Access the first field's value
disp(cMask_path)
disp('is loaded!')
% Eleminate the first segment of all the spokes for accuracies

%
size_cMask = size(cMask);
nbins = size_cMask(1);
cMask = reshape(cMask, [nbins, 22, 2055]); 
cMask(:, 1, :) = []; 

cMask(:, :, 1:p.nShot_off) = []; 
cMask = bmPointReshape(cMask); 


% Run the mitosis function and compute volume elements

[y, t] = bmMitosis(y_tot, t_tot, cMask); 
y = bmPermuteToCol(y); 
ve  = bmVolumeElement(t, 'voronoi_full_radial3' ); 

% Save all the resulting datastructures on the disk. You are now ready
% to run your reconstruction

bmMitosius_create(mDir, y, t, ve); 
disp('Mitosius files are saved!')
disp(mDir)


end

















