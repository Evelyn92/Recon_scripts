% =====================================================
% Author: Yiwei Jia
% Date: Feb 5
% ------------------------------------------------
% initial check the half-resolution result on debi
% Update: concatenate bin with wobin, calculate the heatmap of difference
% =====================================================
%% 1 manually drag the results and then use
bmImage(x_mc);
bmImage(x_wo_mc);

%% 2 roughly confirm roi
% ctrl+shift+y
x_dim = 120/2:180/2;
y_dim = 35:170/2;
z_dim= 90/2:150/2;
x_mc_crop = x_mc{1}(x_dim, y_dim, z_dim);
x_wo_mc_crop = x_wo_mc{1}(x_dim, y_dim, z_dim);
%% 3 check the cropped image
bmImage(x_mc_crop)
%
bmImage(x_wo_mc_crop)
%% Concatenate 2 matrices
x_bin_wobin_concat = cat(2, x_mc_crop, x_wo_mc_crop);
bmImage(x_bin_wobin_concat)

%% diff the 2 matrics
x_diff_crop = x_mc_crop-x_wo_mc_crop;
bmImage(x_diff_crop)
slice_diff = sum(abs(x_diff_crop), [1, 2]); % Sum over rows & columns
[sorted_diff, sorted_indices] = sort(slice_diff(:), 'descend'); % Sort in descending order 
disp(sorted_indices')
% 30    29    33    28    32    27    34    31 / win3 th0.75 cri0.3
% 33    30    46    48    49    32    29    47    31 / win3 th0.98 cri0.3
%% Compare slice_bin and slice_wobin
i = 16;
slice_bin = x_mc_crop(:, :, i);
slice_wobin = x_wo_mc_crop(:, :, i);
compare_two_slices(slice_bin, 'mc', slice_wobin, 'womc', ['Slice-',num2str(i) ])
%% Norm each single slice and see any difference

min_val = 0.00;  % Lower threshold
max_val = 0.98;  % Upper 

slice_bin = x_mc_crop(:, :, i);
slice_bin_norm=norm_slice_mat(slice_bin, min_val, max_val, 'mc', false); %display_flag on

slice_wobin = x_wo_mc_crop(:, :, i);
slice_wobin_norm=norm_slice_mat(slice_wobin, min_val, max_val, 'womc', false);
compare_two_slices(slice_bin_norm, 'mc', slice_wobin_norm, 'womc', ['Norm Slice-',num2str(i) ])
%%

check_intensity_profile(slice_bin_norm,slice_wobin_norm,'mc vs womc');

%%
diff_map_plot(x_mc_crop, x_wo_mc_crop, 16)
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

function diff_map_plot(x1, x2, slice_idx)
x1 = abs(x1);
x2 = abs(x2);
diff_map = abs(x1 - x2);
sgtitle(['Slice', num2str(slice_idx)]);
ax1 = subplot(2, 2, 1);
imagesc(x1(:, :, slice_idx)); axis off; title('mc');
colormap(ax1, gray); % Set colormap only for this subplot
colorbar;
ax2 = subplot(2, 2, 2);
imagesc(x2(:, :, slice_idx)); axis off; title('womc');
colormap(ax2, gray);
colorbar;
ax3 = subplot(2, 2, 3);
imagesc(diff_map(:, :, slice_idx)); axis off; title('Difference');
colormap(ax3, hot); % Only this one uses hot
colorbar;
% ax4 = subplot(2, 2, 4); % Empty subplot for colorbar
% axis off; % Hide axes
% cbar = colorbar(ax3, 'Location', 'westoutside'); % Add colorbar here
% cbar.Label.String = 'Intensity Difference';

end

function check_intensity_profile(img1,img2,info)

imshow(img1, 'InitialMagnification', 'fit'); 
[x, y, intensityProfile1] = improfile;
imshow(img2, 'InitialMagnification', 'fit'); 
[x, y, intensityProfile2] = improfile(img2,x,y);
figure;
plot(intensityProfile1,'-r', 'DisplayName','mc','LineWidth',2);
hold on;
plot(intensityProfile2,'-b', 'DisplayName','womc','LineWidth',2);
legend('show');grid on;
title(['Intensity Profile: ', info]);
end