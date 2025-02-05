clc;
% addpath(genpath("/media/sinf/1,0 TB Disk/Backup/Recon_fork"));
addpath(genpath('/Users/cag/Documents/forclone/Recon_fork'));
base_folder = '/Users/cag/Documents/Dataset/1_Pilot_MREye_Data';
%% Convert the nii to mat.
nii2mat(base_folder, true);
%% Load the .mat files (assuming the variable inside each is named 'dataMatrix')
Sub_num = 4;
cat_num = 5;
Matrix_size = 480;
only_roi = true;
roi_options = {'left', 'right', 'both'};
roi_mode = 'both';

roi_range = cell(1, 4);

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

Bin_title = {'Bin 30177/45210', 'Bin 31466/45210', '', 'Bin 12912/45210' };
woBin_title = {'woBin 42861/45210', 'woBin 42861/45210', '', 'woBin 42861/45210'}; 
woBin_sampled_title = {'woBin sampled 30177/45210','woBin sampled 31466/45210', '', 'woBin sampled 12912/45210'}; 
Discard_title = {'Discard 12831/45210', 'Discard 11556/45210', '', 'Discard 12912/45210'}; 
Vibe_title = 'Vibe';

if Sub_num == 1
    reconDir = '/Users/cag/Documents/Dataset/240922_480/Sub001/';
    %the vibe mat has the opposite direction of libre at z
    %I wanted to align the same slice by
    %finding the first slice from 240-0 direction where the first complete
    %optical nerve occurs -> first_nerve_vibe
    %finding the first slice from 0-240 directio where optical nerve occurs
    %-> first_nerve_libre
    first_nerve_vibe = 128;
    first_nerve_libre = 233;
    slice_offset = first_nerve_libre-(240-first_nerve_vibe);
    % slice_offset = 123;
elseif Sub_num == 2
    reconDir = '/Users/cag/Documents/Dataset/240922_480/Sub002/'; 
    first_nerve_vibe = 125;
    first_nerve_libre = 237;
    slice_offset = first_nerve_libre-(240-first_nerve_vibe);
elseif Sub_num == 3
    reconDir = '/Users/cag/Documents/Dataset/240922_480/Sub003/';
elseif Sub_num == 4
    reconDir = '/Users/cag/Documents/Dataset/240922_480/Sub004/';
    first_nerve_vibe = 136;
    first_nerve_libre = 226;
    slice_offset = first_nerve_libre-(240-first_nerve_vibe);
end

if only_roi
    saveDir = [reconDir, 'cat_6_roi/', roi_mode,'/'];
else
    saveDir = [reconDir, 'cat_6_full/'];
end

if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
%%
for idx = 1:4

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
if cat_num > 4
    base_sub_folder = [base_folder, sprintf('/Sub00%d', Sub_num)];
    if strcmp(roi_mode, roi_options{1})
        % x: vertical from top to the bottom
        % y: horizontal from left to the right
        xEyeRange = 120:240;
        yEyeRange = 1:120;
        roi_range{3} = xEyeRange;
        roi_range{4} = yEyeRange;
    elseif strcmp(roi_mode, roi_options{2})
        % right side
        xEyeRange = 120:240;
        yEyeRange = 120:240;
        roi_range{3} = xEyeRange;
        roi_range{4} = yEyeRange;
    elseif strcmp(roi_mode, roi_options{3})
        % both sides
        xEyeRange = 120:240;
        yEyeRange = 1:240;
        roi_range{3} = xEyeRange;
        roi_range{4} = yEyeRange;
    end

    for idx = 1:cat_num-4
         
         [matName, matDataDir, ~] = uigetfile( ...
                { '*.mat','recon result file (*.mat)'}, ...
                   'Pick a recon result file', ...
                   'MultiSelect', 'off', base_sub_folder);
        
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
          % reverse z axis, and flip the image along x axis
         recon_data = recon_data(:,end:-1:1,end:-1:1);
            
         volume_cell{idx+4} = recon_data;
        
    end

end
%%
for slice_idx = 210:270
% Extract the middle slice along the z-axis for visualization (slice #240 for 480x480x480 matrices)


% Create a figure with 2x2 subplots
figure;

subplot(3, 2, 1);  % First subplot
recon_volume = volume_cell{1};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(Bin_title{Sub_num});
axis equal tight;   % Adjust axis scaling

subplot(3, 2, 2);  % Second subplot
recon_volume = volume_cell{2};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(woBin_title{Sub_num});
axis equal tight;

subplot(3, 2, 3);  % Third subplot
recon_volume = volume_cell{3};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(woBin_sampled_title{Sub_num});
axis equal tight;

subplot(3, 2, 4);  % Fourth subplot
recon_volume = volume_cell{4};
if only_roi
    slice = abs(recon_volume(roi_range{1}, roi_range{2}, slice_idx));
else
    slice = abs(recon_volume(:, :, slice_idx));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);
title(Discard_title{Sub_num});
axis equal tight;


subplot(3, 2, 5);  % Fifth subplot
recon_volume = volume_cell{5};
if only_roi
    slice = abs(recon_volume(roi_range{3}, roi_range{4}, slice_idx-slice_offset));
else
    slice = abs(recon_volume(:, :, slice_idx-slice_offset));
end
slice_norm = mat2gray(double(slice)); % Normalizes the image
imagesc(slice_norm);    % Display first matrix slice
colormap(gray);

title([Vibe_title, sprintf(' Slice-%d', slice_idx-slice_offset)]);
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
