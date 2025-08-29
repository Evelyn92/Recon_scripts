% =====================================================
% Author: Yiwei Jia
% Date: Feb 5
% ------------------------------------------------
% initial check the half-resolution result on debi
% Update: concatenate bin with wobin, calculate the heatmap of difference
% =====================================================
%% 1 manually drag the results and then use
bmImage(x_bin);
bmImage(x_wobin);
% x_bin = flip(x_bin, 1);
% x_wobin = flip(x_wobin, 1);
%% 2 roughly confirm roi
% ctrl+shift+y
x_dim = 150:220;
y_dim = 65:170;
z_dim= 1:240;
x_bin_crop = x_bin(x_dim, y_dim, z_dim);
x_wobin_crop = x_wobin(x_dim, y_dim, z_dim);

% 3 check the cropped image
bmImage(x_bin_crop) 
%
bmImage(x_wobin_crop)
%% Concatenate 2 matrices
x_bin_wobin_concat = cat(2, x_bin_crop, x_wobin_crop);
bmImage(x_bin_wobin_concat)
reconpath = '/Users/cag/Documents/Dataset/recon_results/250423/Sub001/figs_cri_0p08/';
%%

if isfolder(reconpath) == 0
    mkdir(reconpath)
end
savefig = 0;
for idx=100:1:140
    figure;set(gcf, 'Color', 'w');
    subplot(1,1,1);
    imagesc(abs(x_bin_wobin_concat(:, :, idx))); axis off; title(strcat('Bin vs woBin idx: ',num2str(idx))); colorcode = 'gray';colormap(colorcode);   

    ax = gca;
    imgPath = strcat(reconpath, '/bin_vs._woBin-idx_',num2str(idx),'.png');
    if savefig
        disp(['saving: ', imgPath])
        exportgraphics(ax, imgPath)
    end
end
% 
% figure;set(gcf, 'Color', 'w');
% subplot(2,1,1); imagesc(abs(x_bin_wobin_concat(:, :, 225))); axis off; title('Bin vs woBin idx:225'); colorcode = 'gray';colormap(colorcode);   
% subplot(2,1,2); imagesc(abs(x_bin_wobin_concat(:, :, 230))); axis off; title('Bin vs woBin idx:230'); colorcode = 'gray';colormap(colorcode);   
% 
% figure;set(gcf, 'Color', 'w');
% subplot(2,1,1); imagesc(abs(x_bin_wobin_concat(:, :, 260))); axis off; title('Bin vs woBin idx:260'); colorcode = 'gray';colormap(colorcode);   
% subplot(2,1,2); imagesc(abs(x_bin_wobin_concat(:, :, 270))); axis off; title('Bin vs woBin idx:270'); colorcode = 'gray';colormap(colorcode);   
%% diff the 2 matrics
x_diff_crop = x_bin_crop-x_wobin_crop;
bmImage(x_diff_crop)
slice_diff = sum(abs(x_diff_crop), [1, 2]); % Sum over rows & columns
[sorted_diff, sorted_indices] = sort(slice_diff(:), 'descend'); % Sort in descending order 
disp(sorted_indices(1:20)');
%  95    64    94    57    91    58    59    62    87     / win3 th0.75 cri0.3
% 63    64    65    62    60    61    66 / win3 th0.75 cri0.1

%% Compare slice_bin and slice_wobin
i = 66;
slice_bin = x_bin_crop(:, :, i);
slice_wobin = x_wobin_crop(:, :, i);
compare_two_slices(slice_bin, 'bin', slice_wobin, 'wobin', ['Slice-',num2str(i) ])
% Norm each single slice and see any difference

min_val = 0.00;  % Lower threshold
max_val = 0.8;  % Upper threshold

slice_bin = x_bin_crop(:, :, i);
slice_bin_norm=norm_slice_mat(slice_bin, min_val, max_val, 'bin', false); %display_flag on

slice_wobin = x_wobin_crop(:, :, i);
slice_wobin_norm=norm_slice_mat(slice_wobin, min_val, max_val, 'wobin', false);
compare_two_slices(slice_bin_norm, 'bin', slice_wobin_norm, 'wobin', ['Norm Slice-',num2str(i) ])

%%
function compare_two_slices(slice1, slice1_label, slice2, slice2_label, sg_label)
    slice1 = mat2gray(abs(slice1));
    slice2 = mat2gray(abs(slice2));
    figure;
 
    h = imshow(cat(2, slice1, slice2), 'InitialMagnification', 'fit'); 
    ax = gca; % Get the current axes
    axis tight;
    set(ax, 'FontName', 'Times New Roman', 'FontSize', 20);  % Set font for axis
    title(sg_label, 'FontSize', 20, 'FontName', 'Times New Roman');
    % title(sg_label,'FontSize',20, 'FontName','Times New Roman');
    colorcode = 'viridis(256)';
    colormap(colorcode);   
    colorbar;          % Show color scale
end

function slice_redistributed = norm_slice_mat(org_slice, min_val, max_val, slice_label, display_flag)
slice_norm = mat2gray(abs(org_slice));
% Clip the intensity values
slice_I_clipped = slice_norm;
slice_I_clipped(slice_norm < min_val) = min_val;
slice_I_clipped(slice_norm > max_val) = max_val;

% Normalize to [0, 1] after clipping
slice_redistributed = (slice_I_clipped - min(slice_I_clipped(:))) / (max(slice_I_clipped(:)) - min(slice_I_clipped(:)));

if display_flag
    % Display the results
    figure;sgtitle(slice_label);
    subplot(1,2,1); imshow(slice_norm, []); title('Original Image');
    subplot(1,2,2); imshow(slice_redistributed, []); title('Clipped & Redistributed Image');
    

    % figure;
    % subplot(1,2,1); imhist(slice_norm); title('Original Histogram');
    % subplot(1,2,2); imhist(slice_redistributed); title('Clipped & Redistributed Histogram');
end

end
%%
close all;
for slice_idx=100:140
    save_fig=1;
    diff_map_plot(x_bin_crop, x_wobin_crop, slice_idx, 0, reconpath)
end
function diff_map_plot(x1, x2, slice_idx, save_fig, reconpath)
addpath(genpath('/Users/cag/Documents/forclone/mfuderer-colorResources-05c224c'));
x1 = abs(x1);
x2 = abs(x2);
diff_map = abs(x1 - x2);
figure; set(gcf, 'color', 'w')

ax1 = subplot(3, 1, 1);
imshow(x1(:, :, slice_idx), [0,3]); axis off; title(strcat('bin-', num2str(slice_idx)));
colormap(ax1, gray); % Set colormap only for this subplot
colorbar;

ax2 = subplot(3, 1, 2);
% imagesc
imshow(x2(:, :, slice_idx), [0,3]); axis off;  title(strcat('wobin-', num2str(slice_idx)));
colormap(ax2, gray);
colorbar;

ax3 = subplot(3, 1, 3);
loLev = 0.05;upLev = 0.8;[imClip, rgb_vec] = relaxationColorMap('T1',diff_map(:,:,slice_idx) , loLev, upLev);
imshow(imClip, 'DisplayRange', [loLev, upLev], 'InitialMagnification', 'fit'); title(strcat('difference map-', num2str(slice_idx)));
colormap(ax3, rgb_vec); colorbar;


if isfolder(reconpath) == 0
    mkdir(reconpath)
end
if save_fig
imgPath = strcat(reconpath, '/diffMap-idx_',num2str(slice_idx),'.png');
exportgraphics(gcf, imgPath, 'Resolution', 300);  % high-res PNG
end
% imagesc(diff_map(:, :, slice_idx)); axis off; title('Difference');
% colormap(ax3, hot); % Only this one uses hot
% colorbar;
% ax4 = subplot(2, 2, 4); % Empty subplot for colorbar
% axis off; % Hide axes
% cbar = colorbar(ax3, 'Location', 'westoutside'); % Add colorbar here
% cbar.Label.String = 'Intensity Difference';

end