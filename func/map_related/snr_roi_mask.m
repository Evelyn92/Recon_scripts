close all;clc;
%% please load your recon image here
% xrms is a 3D volume
sliceIdx = 120; % select the slice index
imgSlice_idea =xrms_idea(:,:,sliceIdx);
imgSlice_pulseq =xrms_pulseq(:,:,sliceIdx);

imgSlice_idea = normImage(imgSlice_idea);
imgSlice_pulseq = normImage(imgSlice_pulseq);
%%
SelectBack=1;
CalSNR = 1; % 1: calculate SNR -- 0: no SNR


cc = 1; % set the phantom circle index
imgIdx = 1; % set the image index for defining note suffix
note_suffix = strcat('cc', num2str(cc), '_sl', num2str(sliceIdx), '_img', num2str(imgIdx)); % cc1 for 'circle 1'
roiShape = 'rect'; %option: 'circle' 'rect'
cal_roi_mask(imgSlice, roiShape, SelectBack, CalSNR, note_suffix);
%% You can check the overlay here

check_overlay_mask_img(maskBg, imgSlice_idea);
check_overlay_mask_img(maskBg, imgSlice_pulseq);
%%
close all;
% reorder the mask from 1 to 18
maskROI_array = {maskROI_cc15_sl120_img1, maskROI_cc11_sl120_img1, ...
    maskROI_cc5_sl120_img1, maskROI_cc3_sl120_img1, ...
    maskROI_cc17_sl120_img1, maskROI_cc13_sl120_img1, ...
    maskROI_cc7_sl120_img1, maskROI_cc1_sl120_img1,...
    maskROI_cc18_sl120_img1, maskROI_cc14_sl120_img1, ...
    maskROI_cc9_sl120_img1, maskROI_cc6_sl120_img1,...
    maskROI_cc2_sl120_img1, maskROI_cc16_sl120_img1,...
    maskROI_cc12_sl120_img1, maskROI_cc10_sl120_img1,...
    maskROI_cc8_sl120_img1, maskROI_cc4_sl120_img1};

snr_estimate_list_idea = {};
for maskIndex = 1: length(maskROI_array) 
    maskROI = maskROI_array{maskIndex};
    
    check_overlay_mask_img(maskROI, imgSlice_idea);
    roi_values = imgSlice_idea(maskROI);
    bg_values = imgSlice_idea(maskBg);
    snr_estimate = calSNR(roi_values, bg_values);
    snr_estimate_list_idea{maskIndex} = snr_estimate;
end

clear roi_values bg_values snr_estimate;
snr_estimate_list_pulseq = {};
for maskIndex = 1: length(maskROI_array) 
    maskROI = maskROI_array{maskIndex};
    
    check_overlay_mask_img(maskROI, imgSlice_pulseq);

    roi_values = imgSlice_pulseq(maskROI);
    bg_values = imgSlice_pulseq(maskBg);
    snr_estimate = calSNR(roi_values, bg_values);
    snr_estimate_list_pulseq{maskIndex} = snr_estimate;
end
clear roi_values bg_values snr_estimate;
%%
snr1 = cell2mat(snr_estimate_list_idea);
snr2 = cell2mat(snr_estimate_list_pulseq);

%% Paired line plot
colors = lines(18);  % or use 'parula(18)', 'jet(18)', etc.
figure; set(gcf, 'Color', 'w');
hold on;
for i = 1:18
    plot([1 2], [snr1(i), snr2(i)], '-o', 'Color', colors(i,:), 'LineWidth', 1.5);
end
xticks([1 2]);
xticklabels({'IDEA LIBRE 6p2', 'pulseq LIBRE 6p2'});
ylabel('SNR');
title('Paired SNR Comparison');
grid on;

% Optional: show group means
plot([1 2], [mean(snr1), mean(snr2)], 'k*-', 'LineWidth', 2);
legend('Individual pairs','Group mean');
%%
% Plot
figure; hold on;
plot(1:18, snr1, '-o', 'LineWidth', 2, 'DisplayName', 'IDEA LIBRE 6p2' );
plot(1:18, snr2, '-s', 'LineWidth', 2, 'DisplayName', 'pulseq LIBRE 6p2');

% Beautify
xlabel('Subject / Pair Index');
ylabel('SNR');
title('SNR Comparison Across Pairs');
legend('Location', 'best');
grid on;
set(gcf, 'Color', 'w');  % white figure background
%%
function cal_roi_mask(imageSlice, roiShape, SelectBack, CalSNR, note_suffix)
    % imageSlice: 2D matrix (e.g., one MRI slice)
   if CalSNR
       SelectBack = 1;
       warning('SelectBack should be 1 for SNR calculation!')
   end
   fprintf('select ROI!')
   [roi_values, maskROI] = select_ROI(imageSlice, roiShape);
   maskName =  strcat('maskROI_', note_suffix);

   assignin('base',maskName, maskROI);  % Optional: export to workspace
   assignin('base', 'roi_values', roi_values);  % Optional: export to workspace
   if SelectBack
       fprintf('select background!')
       [bg_values, maskBg] = select_ROI(imageSlice, 'rect');
       assignin('base', 'maskBg', maskBg);  % Optional: export to workspace
       fprintf('The background mask maskBg has been saved to workspace! \n')
       assignin('base', 'bg_values', bg_values);  % Optional: export to workspace
       fprintf('The background mean value bg_values has been saved to workspace! \n')
   end
   
   if SelectBack && CalSNR
      snr_estimate = calSNR(roi_values, bg_values);
      assignin('base', 'snr_estimate', snr_estimate);  % Optional: export to workspace
      fprintf('The SNR snr_estimate %f has been saved to workspace!', snr_estimate);
   end 
    
end


function [roi_values, mask]= select_ROI(imageSlice, shapeInfo)
    figure;
    imshow(imageSlice, [], 'InitialMagnification', 'fit');

    if strcmp(shapeInfo, 'circle')
        title('Draw a circle ROI. Double-click to confirm.'); 
        % Let user select circle ROI
        h = imellipse(gca, []);
    elseif strcmp(shapeInfo, 'rect')
        title('Draw a rect ROI. Double-click to confirm.'); 
        % Let user select rectangular background
        h = imrect(gca, []);
       
    end

    % Create binary mask from the selected region
    position = wait(h);  % Wait for user to finish drawing
    mask = createMask(h); 

    check_overlay_mask_img(mask, imageSlice);
    roi_values = imageSlice(mask);

end

function check_overlay_mask_img(mask, imageSlice)
    % Show the mask overlay
    figure;
    imshow(imageSlice, [], 'InitialMagnification', 'fit'); hold on;
    visboundaries(mask, 'Color', 'r');
    title('Selected ROI Overlay');

end

function snr_estimate = calSNR(roi_values, bg_values)
        snr_estimate = mean(roi_values) / std(bg_values);
       fprintf('Estimated SNR in ROI: %.2f\n', snr_estimate);
end

function normedSlice = normImage(imgSlice)
    mag = abs(imgSlice);
    normedSlice = (mag - min(mag(:)))/(max(mag(:)) - min(mag(:)));
    disp('Image is normalized to [0,1]')
end