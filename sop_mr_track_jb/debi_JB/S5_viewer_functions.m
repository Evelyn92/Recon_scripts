clear all;
clc;

addpath(genpath('/media/debi/1,0 TB Disk/Backup/Recon_fork'));
%% Load the .mat files (assuming the variable inside each is named 'dataMatrix')
Sub_num = 1;
cat_num = 4;
Matrix_size = 240;
only_roi = true;
roi_options = {'left', 'right', 'both'};
roi_mode = 'both';

roi_range = cell(1, 2);

if strcmp(roi_mode, roi_options{1})
%left side
xLeftEyeRange = (Matrix_size-Matrix_size/2):(Matrix_size-Matrix_size/4);
yLeftEyeRange = (Matrix_size-Matrix_size*3/4):(Matrix_size-Matrix_size/2);
roi_range{1} = xLeftEyeRange;
roi_range{2} = yLeftEyeRange;
elseif strcmp(roi_mode, roi_options{2})
%right side
xRightEyeRange = (Matrix_size-Matrix_size/2):(Matrix_size-Matrix_size/4);
yRightEyeRange = (Matrix_size-Matrix_size/2):(Matrix_size-Matrix_size/4);
roi_range{1} = xRightEyeRange;
roi_range{2} = yRightEyeRange;
elseif strcmp(roi_mode, roi_options{3})
%zoom in the eye region
xBothEyesRange = (Matrix_size-Matrix_size/2):(Matrix_size-Matrix_size/4);
yBothEyesRange = (Matrix_size-Matrix_size*3/4):(Matrix_size-Matrix_size/4);
roi_range{1} = xBothEyesRange;
roi_range{2} = yBothEyesRange;
end

volume_cell = cell(cat_num, 1);

% Bin_title = {'Bin 30177/45210', 'Bin 31466/45210', '', 'Bin 12912/45210' };
% woBin_title = {'woBin 42861/45210', 'woBin 42861/45210', '', 'woBin 42861/45210'}; 
% woBin_sampled_title = {'woBin sampled 30177/45210','woBin sampled 31466/45210', '', 'woBin sampled 12912/45210'}; 
% Discard_title = {'Discard 12831/45210', 'Discard 11556/45210', '', 'Discard 12912/45210'}; 
reconDir = '/media/debi/1,0 TB Disk/241120_JB/Sub001/';
plot_title={'woBin', 'left', 'right', 'up'};

if only_roi
    saveDir = [reconDir, 'cat_4_roi/', roi_mode,'/'];
else
    saveDir = [reconDir, 'cat_4_full/'];
end

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
%%
for idx = 1:cat_num

 [matName, matDataDir, ~] = uigetfile( ...
        { '*.mat','recon result file (*.mat)'}, ...
           'Pick a recon result file', ...
           'MultiSelect', 'off', reconDir);

    if matName == 0
        warning('No mask file selected');
        return;
    end
    filepathMatData = fullfile(matDataDir, matName);
    disp(filepathMatData);

    recon_data_cell = load(filepathMatData);
    fieldNames = fieldnames(recon_data_cell); 
    recon_data_cell = recon_data_cell.(fieldNames{1});
    if iscell(recon_data_cell)
        recon_data = recon_data_cell{1};
    else
        recon_data = recon_data_cell;
    end


 volume_cell{idx} = recon_data;

end

% Perform pairwise comparison
volume1 = volume_cell{1};
volume2 = volume_cell{2};
volume3 = volume_cell{3};
volume4 = volume_cell{4};
isEqual12 = isequal(volume1, volume2);
isEqual13 = isequal(volume1, volume3);
isEqual14 = isequal(volume1, volume4);
isEqual23 = isequal(volume2, volume3);
isEqual24 = isequal(volume2, volume4);
isEqual34 = isequal(volume3, volume4);

% Check if any of the comparisons resulted in equality
if isEqual12 || isEqual13 || isEqual14 || isEqual23 || isEqual24 || isEqual34
    disp('Some volumes are equal. Check the volumes.');
else
    disp('All volumes are different.');
end

%%
for slice_idx = 100:150
% Extract the middle slice along the z-axis for visualization (slice #240 for 480x480x480 matrices)


% Create a figure with 2x2 subplots
figure;

subplot(2, 2, 1);  % First subplot
recon_volume = volume_cell{1};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(plot_title{1});
axis equal tight;   % Adjust axis scaling

subplot(2, 2, 2);  % Second subplot
recon_volume = volume_cell{2};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(plot_title{2});
axis equal tight;

subplot(2, 2, 3);  % Third subplot
recon_volume = volume_cell{3};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(plot_title{3});
axis equal tight;

subplot(2, 2, 4);  % Fourth subplot
recon_volume = volume_cell{4};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(plot_title{4});
axis equal tight;

sgtitle(sprintf('Subject %d - Slice %d', Sub_num, slice_idx));
saveas(gcf, [saveDir, sprintf('Slice_%d.png', slice_idx)]);
disp(['All the slices are saved here', saveDir]);
end
% % Add a colorbar to each subplot (optional)
% for i = 1:4
%     subplot(2, 2, i);
%     colorbar;
% end
