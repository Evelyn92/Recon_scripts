function concatenate_slices(Binning_png, woBinning_png, concat_output_dir, slice_start, slice_end, only_roi,roi_range)
    %% concatenation and inspect
    % concat_output_dir = '/media/sinf/1,0 TB Disk/240829_recon_for_poster/Sub002/concat';

    if ~exist(concat_output_dir, 'dir')
            % If it doesn't exist, create it
            mkdir(concat_output_dir);
            disp(['Directory created: ', concat_output_dir]);
    else
        disp(['Directory already exists: ', concat_output_dir]);
    end
    if only_roi
        xRange = roi_range{1};
        yRange = roi_range{2};
        disp('xRange')
        disp([num2str(xRange(1)),' ', num2str(xRange(end))])
        disp('yRange')
        disp([num2str(yRange(1)),' ', num2str(yRange(end))])
    end

    for k = slice_start:slice_end
        img_1 = imread(fullfile(Binning_png, sprintf('slice_%03d.png', k)));
        img_2 = imread(fullfile(woBinning_png, sprintf('slice_%03d.png', k)));
        % img_3 = imread(fullfile(sub002_t1vibe_png, sprintf('slice_%03d.png', k)));
        % img_3_rot = imrotate(img_3, -90);
       
        % imshow(img_k_scale(112:172,20:70));
        % up:down, left:right
        if only_roi
            concatenated_image = [img_1(xRange, yRange), img_2(xRange, yRange)];%%, img_3_rot(60:120,8:58) up:down, left:right
        else
            concatenated_image = [img_1, img_2];%% , img_3] up:down, left:right
        end
        % Display the concatenated image
        imshow(concatenated_image);
        % Optionally save the concatenated image
        concate_img_file = fullfile(concat_output_dir, sprintf('slice_%03d.png', k));
        imwrite(concatenated_image, concate_img_file);
    end
    
    disp('Finish the concatenation!')
    disp(concat_output_dir)
end