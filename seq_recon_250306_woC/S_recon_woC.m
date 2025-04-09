clear;clc; 
% =====================================================
% Author: Yiwei Jia
% Date: March 10
% ------------------------------------------------
% Recon without C with Mathilda
% Update: just for test before running on HPC


%% =====================================================

subject_num = 3;
datasetDir = '/home/debi/yiwei/mreye_dataset/250314_libre_pulseq/';
raw_data = fullfile(datasetDir,'meas_MID00150_FID291941_d_t1w_rect_rf.dat');
reconDir = '/home/debi/yiwei/recon_results/250314_libre/';

mDir = [reconDir, '/Sub003/T1_LIBRE_woBinning/mitosius/mask_', 'noC_rect', '/'];

%
y   = bmMitosius_load(mDir, 'y'); 
t   = bmMitosius_load(mDir, 't'); 
ve  = bmMitosius_load(mDir, 've'); 

disp('Mitosius has been loaded!')
%% compileScript()
Matrix_size = 240;
ReconFov = 480; %mm
N_u     = [Matrix_size, Matrix_size, Matrix_size]; % Matrix size: Size of the Virtual cartesian grid in the fourier space (regridding)
n_u     = [Matrix_size, Matrix_size, Matrix_size]; % Image size (output)
dK_u    = [1, 1, 1]./ReconFov; % Spacing of the virtual cartesian grid
nFr     = 1; 
% best achivable resolution is 1/ N_u*dK_u If you have enough coverage
%%
nCh = size(y{1}, 2);
x0 = cell(nCh, 1);
for i = 1:nFr
    for iCh = 1:nCh
    x0{iCh} = bmMathilda(y{i}(:,iCh), t{i}, ve{i}, [], N_u, n_u, dK_u, [], [], [], []);
    disp(['Processing channel: ', num2str(iCh),'/', num2str(nCh)])
    end
end

%
bmImage(x0);
%%
mask_note = 'noC_rect';
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
