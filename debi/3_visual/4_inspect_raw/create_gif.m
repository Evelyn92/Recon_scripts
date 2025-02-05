%% Data loader
clc;
%%
mode = 'raw';
% 'lc' 'raw''clean'

subject_num = 1;
delta = 1;

%%
if strcmp(mode, 'clean')
    x={};
if subject_num == 1
    if delta == 0.1
        x_ = load(['/Users/cag/Documents/Dataset/241120_JB/' ...
        'Sub001/T1_LIBRE_Binning/output/mask_0/x_nIter20_delta_0.10.mat']);
        x{1,1} = x_.x;
        
        x_ = load(['/Users/cag/Documents/Dataset/241120_JB/' ...
            'Sub001/T1_LIBRE_Binning/output/mask_1/x_nIter20_delta_0.10.mat']);
        x{2,1}=x_.x;
        
        x_ =load(['/Users/cag/Documents/Dataset/241120_JB/' ...
            'Sub001/T1_LIBRE_Binning/output/mask_2/x_nIter20_delta_0.10.mat']);
        x{3,1}=x_.x;
        
        x_ =load(['/Users/cag/Documents/Dataset/241120_JB/' ...
            'Sub001/T1_LIBRE_Binning/output/mask_3/x_nIter20_delta_0.10.mat']);
        x{4,1}=x_.x;
    else
        disp('aaa')
        x_ = load(['/Users/cag/Documents/Dataset/241120_JB/' ...
        'Sub001/T1_LIBRE_Binning/output/mask_0/x_nIter20_delta_1.00.mat']);
        x{1,1} = x_.x;
        
        x_ = load(['/Users/cag/Documents/Dataset/241120_JB/' ...
            'Sub001/T1_LIBRE_Binning/output/mask_1/x_nIter20_delta_1.00.mat']);
        x{2,1}=x_.x;
        
        x_ =load(['/Users/cag/Documents/Dataset/241120_JB/' ...
            'Sub001/T1_LIBRE_Binning/output/mask_2/x_nIter20_delta_1.00.mat']);
        x{3,1}=x_.x;
        
        x_ =load(['/Users/cag/Documents/Dataset/241120_JB/' ...
            'Sub001/T1_LIBRE_Binning/output/mask_3/x_nIter20_delta_1.00.mat']);
        x{4,1}=x_.x;

    end

elseif subject_num == 2
else %subject_num == 3
    if delta == 0.1
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_0/x_nIter20_delta_0.10.mat']);
        x{1,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_1/x_nIter20_delta_0.10.mat']);
        x{2,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_2/x_nIter20_delta_0.10.mat']);
        x{3,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_3/x_nIter20_delta_0.10.mat']);
        x{4,1} = x_.x;
    else
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_0/x_nIter20_delta_1.00.mat']);
        x{1,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_1/x_nIter20_delta_1.00.mat']);
        x{2,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_2/x_nIter20_delta_1.00.mat']);
        x{3,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/mask_3/x_nIter20_delta_1.00.mat']);
        x{4,1} = x_.x;
    end

end

end
x_clean = x;
%%
if strcmp(mode, 'raw')
   
x={};
if subject_num == 1
    if delta==0.1
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_0/x_nIter20_delta_0.10.mat');
        x{1,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_1/x_nIter20_delta_0.10.mat');
        x{2,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_2/x_nIter20_delta_0.10.mat');
        x{3,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_3/x_nIter20_delta_0.10.mat');
        x{4,1} = x_.x;
    else
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_0/x_nIter20_delta_1.00.mat');
        x{1,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_1/x_nIter20_delta_1.00.mat');
        x{2,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_2/x_nIter20_delta_1.00.mat');
        x{3,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/raw_mask_3/x_nIter20_delta_1.00.mat');
        x{4,1} = x_.x;
    end

elseif subject_num == 3
    if delta == 0.1
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_0/x_nIter20_delta_0.10.mat']);
        x{1,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_1/x_nIter20_delta_0.10.mat']);
        x{2,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_2/x_nIter20_delta_0.10.mat']);
        x{3,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_3/x_nIter20_delta_0.10.mat']);
        x{4,1} = x_.x;
    else
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_0/x_nIter20_delta_1.00.mat']);
        x{1,1} = x_.x;

        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_1/x_nIter20_delta_1.00.mat']);
        x{2,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_2/x_nIter20_delta_1.00.mat']);
        x{3,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/raw_mask_3/x_nIter20_delta_1.00.mat']);
        x{4,1} = x_.x;
    end

end
end
x_raw = x;

%%
if strcmp(mode, 'lc')
   
x={};
if subject_num == 1
    if delta==0.1
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_0/x_nIter20_delta_0.10.mat');
        x{1,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_1/x_nIter20_delta_0.10.mat');
        x{2,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_2/x_nIter20_delta_0.10.mat');
        x{3,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_3/x_nIter20_delta_0.10.mat');
        x{4,1} = x_.x;
    else
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_0/x_nIter20_delta_1.00.mat');
        x{1,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_1/x_nIter20_delta_1.00.mat');
        x{2,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_2/x_nIter20_delta_1.00.mat');
        x{3,1} = x_.x;
        x_ = load('/Users/cag/Documents/Dataset/241120_JB/Sub001/T1_LIBRE_Binning/output/lc_mask_3/x_nIter20_delta_1.00.mat');
        x{4,1} = x_.x;
    end

elseif subject_num == 3
    if delta == 0.1
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_0/x_nIter20_delta_0.10.mat']);
        x{1,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_1/x_nIter20_delta_0.10.mat']);
        x{2,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_2/x_nIter20_delta_0.10.mat']);
        x{3,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_3/x_nIter20_delta_0.10.mat']);
        x{4,1} = x_.x;
    else
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_0/x_nIter20_delta_1.00.mat']);
        x{1,1} = x_.x;

        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_1/x_nIter20_delta_1.00.mat']);
        x{2,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_2/x_nIter20_delta_1.00.mat']);
        x{3,1} = x_.x;
    
        x_ = load(['/Volumes/MatTechLab/Yiwei_recon/241120_JB/' ...
            'Sub003/T1_LIBRE_Binning/output/lc_mask_3/x_nIter20_delta_1.00.mat']);
        x{4,1} = x_.x;
    end

end
end
x_lc=x;
%% equalize 3d matrix
%remember change the name
x_eq = cell(4,1);
for fr=1:4
    img3D = abs(x{fr,1});
    equalized3D=equalize_image_3d(abs(img3D));
    x_eq{fr,1} = equalized3D;
end

%%
% data = load('your_file.mat');
% Assume 'frames' is the variable in the .mat file containing the image frames
frames = data.frames;

gif_filename = 'output.gif'; % Name of the output GIF file
delay_time = 0.1; % Time delay between frames in seconds

for i = 1:size(frames, 3) % Assuming frames is a 3D matrix (height x width x frames)
    % Extract the current frame
    current_frame = frames(:, :, i);
    
    % Normalize the frame if needed (e.g., grayscale images should be scaled 0-255)
    if ~isinteger(current_frame)
        current_frame = uint8(255 * mat2gray(current_frame)); % Convert to uint8
    end
    
    % Write to the GIF
    if i == 1
        % First frame, create the GIF file
        imwrite(current_frame, gif_filename, 'gif', 'LoopCount', Inf, 'DelayTime', delay_time);
    else
        % Subsequent frames, append to the GIF
        imwrite(current_frame, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

function equalized3D=equalize_image_3d(img3D)
    % Assume img3D is your 3D image matrix (e.g., size: [rows, cols, slices])
    [rows, cols, slices] = size(img3D);
    equalized3D = zeros(size(img3D)); % Initialize the output matrix
    
    for i = 1:slices
        % Extract the ith 2D slice
        slice = img3D(:, :, i);
        slice_norm = mat2gray(slice); % Normalizes the image
        min_val = 0.01;  % Lower threshold
        max_val = 0.8;  % Upper threshold
        slice_I_clipped = slice_norm;
        slice_I_clipped(slice_norm < min_val) = min_val;
        slice_I_clipped(slice_norm > max_val) = max_val;
        slice_redistributed = (slice_I_clipped - min(slice_I_clipped(:))) / (max(slice_I_clipped(:)) - min(slice_I_clipped(:)));
        equalizedSlice = slice_redistributed;
        
        % Store the equalized slice back
        equalized3D(:, :, i) = equalizedSlice;
    end
    
    % equalized3D now contains the histogram-equalized 3D image
end