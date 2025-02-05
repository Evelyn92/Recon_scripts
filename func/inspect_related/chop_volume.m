function chop_volume(input_file, output_dir, only_roi, roi_range, equal)
    % Create the output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        % If it doesn't exist, create it
        mkdir(output_dir);
        disp(['Directory created: ', output_dir]);
    else
        disp(['Directory already exists: ', output_dir]);
    end

    if only_roi
        xRange = roi_range{1};
        yRange = roi_range{2};
        disp('xRange')
        disp([num2str(xRange(1)),' ', num2str(xRange(end))])
        disp('yRange')
        disp([num2str(yRange(1)),' ', num2str(yRange(end))])
    end
        
    % Determine file extension
    [~, ~, ext] = fileparts(input_file);
    
    % Check the file extension and process accordingly
    switch ext
        case '.nii'
            disp('processing nifti files')
            % Load the NIfTI file
            recon_data = niftiread(input_file);

            % Get the size of the 3D volume
            
            [~, ~, slices] = size(recon_data);
           
            
            % Loop through each slice along the z-axis (axial slices)
            for i = 1:slices
                % Extract the i-th slice
                if only_roi
                    slice = abs(recon_data(:, :, i));
                    rotated_slice = rot90(slice, -1);% -1 means 90 degrees clockwise
                    slice = abs(rotated_slice(xRange, yRange));
                else
                    slice = abs(recon_data(:, :, i));
                    rotated_slice = rot90(slice, -1);% -1 means 90 degrees clockwise
                    slice = abs(rotated_slice);
                end
                
                % Normalize the slice to the range [0, 1] for saving as image
                slice_norm = mat2gray(slice); % Normalizes the image
                % equalize
                if equal
                    min_val = 0.01;
                    max_val = 0.8;
                    slice_equal = equal_func(slice_norm, min_val, max_val);
                else 
                    slice_equal = slice_norm;
                end
              
                % Create a filename for each slice
                filename = fullfile(output_dir, sprintf('slice_%03d.png', i));
                imshow(slice_equal);
                % Save the slice as an image
                imwrite(slice_equal, filename);
            end
            
            disp('All slices have been saved from nii.');

        otherwise
            disp('processing .mat files')
            recon_data_cell = load(input_file);
            fieldNames = fieldnames(recon_data_cell); 
            recon_data_cell = recon_data_cell.(fieldNames{1});
            if iscell(recon_data_cell)
                recon_data = recon_data_cell{1};
            else
                recon_data = recon_data_cell;
            end
            % Get the size of the 3D volume
            [~, ~, slices] = size(recon_data);
            % Loop through each slice along the z-axis (axial slices)
            for i = 1:slices
                % Extract the i-th slice
                if only_roi
                    slice = abs(recon_data(xRange, yRange, i));
                else
                    slice = abs(recon_data(:, :, i));
                end

                
                % Normalize the slice to the range [0, 1] for saving as image
                slice_norm = mat2gray(double(slice)); % Normalizes the image
                % % equalize
                if equal
                    slice_equal = histeq(slice_norm);
                else
                    slice_equal = slice_norm;
                end
         
                % Create a filename for each slice
                filename = fullfile(output_dir, sprintf('slice_%03d.png', i));
                
                % Save the slice as an image
                imwrite(slice_equal, filename);
               
            end
            disp('All slices have been saved from mat.');
            disp(output_dir);
      
    end
    
    



end

