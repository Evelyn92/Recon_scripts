
addpath(genpath('/Users/cag/Documents/forclone/Recon_scripts'));
addpath(genpath('/Users/cag/Documents/forclone/pulseq'));
addpath(genpath('/Users/cag/Documents/forclone/monalisa'));
%play with 0822 sub002
subject_num = 3;


datasetDir = '/Users/cag/Documents/Dataset/datasets/250822/';
reconDir = '/Users/cag/Documents/Dataset/recon_results/250822/playRecon/';
if subject_num == 2
    % gx = -Gx, gy=-Gy, gz=-Gz
    % meas_MID00158_FID314157_t1w_swap_0.dat
    meas_name_suffix = '_MID00158_FID314157_t1w_swap_0';
    mask_note = 'full_radial3_phylotaxis';%'full_radial3_phylotaxis_flipxyz';%'full_radial3_phylotaxis';
elseif subject_num == 3
    meas_name_suffix = '_MID00159_FID314158_t1w_swap_1';
    mask_note = 'full_radial3_phylotaxis_swapxy';%'full_radial3_phylotaxis_swapxy';%'full_radial3_phylotaxis'
end

nShot=1000;
meas_name = ['meas', meas_name_suffix];
measureFile = [datasetDir, meas_name,'.dat'];
reader = createRawDataReader(measureFile, true);
reader.acquisitionParams.nSeg = 22;
reader.acquisitionParams.nShot = 1000;
reader.acquisitionParams.nShot_off = 14;


% play with this part

reader.acquisitionParams.traj_type = mask_note;

% Ensure consistency in number o1f shot-off points
nShotOff = reader.acquisitionParams.nShot_off;
% Acquisition from Bern need to change the following part!!
p = reader.acquisitionParams;
p.nShot_off = 14; % in case no validation UI
p.nShot = nShot; % in case no validation UI
p.nSeg = 22; % in case no validation UI
% Load the raw data and compute trajectory and volume elements
y_tot = reader.readRawData(true, true);  % Filter nshotoff and SI
t_tot = bmTraj(p);                       % Compute trajectory
ve_tot = bmVolumeElement(t_tot, 'voronoi_full_radial3');  % Volume elements
% Adjust grid size for coil sensitivity maps
FoV = p.FoV;  % Field of View

voxel_size = 2; % so matrix_size = 120 now
matrix_size = round(FoV/voxel_size);  % Max nominal spatial resolution
N_u = [matrix_size, matrix_size, matrix_size];
dK_u = [1, 1, 1]./FoV;

nCh = size(y_tot, 1);
nCh
nFr = 1;
x0 = cell(nCh, 1);
for i = 1:nFr
    for iCh = 1:nCh
    x0{iCh} = bmMathilda(y_tot(iCh,:), t_tot, ve_tot, [], N_u, N_u, dK_u, [], [], [], []);
    disp(['Processing channel: ', num2str(iCh),'/', num2str(nCh)])
   
    end
end
bmImage(x0);
x0Dir = [reconDir, '/Sub00',num2str(subject_num),'/traj_',mask_note,'/'];
 
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
% Root mean square across the channels
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
bmImage(xrms)