clc; clear;
addpath(genpath("/Users/cag/Documents/forclone/monalisa/"))
addpath(genpath("/Users/cag/Documents/forclone/mapVBVD_Jaime/"))
%
%% ==========================================================
% Author: Yiwei Jia
% Date: May 2, 2025
% Here comes brain and binning!
% derived from /Users/cag/Documents/forclone/mapVBVD/twix_process_yj
% --------------------------------------------------------------------------------------------------
% We need to intially inspect MR raw data for
% (1) QC of raw data once a new batch of acquisition is done
% (2) Check the duration of raw data to determine the length of ET mask
% (3) Determine the start of trigger from physio: Now scanner -> trigger
% ============================================================

subject_num = 4;
if subject_num == 6
    meas_name_suffix = '_MID00022_FID298741_JB_LIBRE2p2_a8_woPERewinder';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;

elseif subject_num == 7
    meas_name_suffix = '_MID00023_FID298742_JB_RectPulse_a5_woPERewinder';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;

elseif subject_num == 1
    meas_name_suffix = '_MID00031_FID298750_sub1_debug';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;

elseif subject_num == 2
    meas_name_suffix = '_MID00032_FID298751_sub2_debug';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;
elseif subject_num == 3
    meas_name_suffix = '_MID00033_FID298752_sub3_debug';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;

elseif subject_num == 4
    meas_name_suffix = '_MID00034_FID298753_sub4_debug';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;

elseif subject_num == 5
    meas_name_suffix = '_MID00035_FID298754_sub5_debug';
    hc_name_suffix = '_MID00030_FID298749_sub1_HC';
    bc_name_suffix = '_MID00028_FID298747_sub1_BC';
    nShot=1000;
end

% meas_name = ['meas', meas_name_suffix];
meas_name = ['meas', meas_name_suffix];
hc_name = ['meas', hc_name_suffix];
bc_name = ['meas', bc_name_suffix];

twix_name = ['twix', meas_name_suffix];
twix_path =  ['/Users/cag/Documents/Dataset/datasets/250502/', twix_name,'.mat'];
datasetDir = '/Users/cag/Documents/Dataset/datasets/250502/';
raw_data = [datasetDir, meas_name,'.dat'];
hc_data = [datasetDir, hc_name,'.dat'];
bc_data = [datasetDir, bc_name,'.dat'];
save_trigger = false;
% (1) (2)
inspect_raw_data = true;
if inspect_raw_data
    bmTwix_info(raw_data);
    % check_data_timestamp(datasetDir,false);
end
% Please paste the output here for better inspection
%-----------------------------------------------------------------
% In protocol, we set the T1 LIBRE duration: 362(sec)*1000(freq)
% And we set the T1 VIBE duration: 357(sec)*1000(freq)
% However, here the duration of the rawdata is: 361672.5 ms with data points:45210
%-----------------------------------------------------------------

%% (3)
twix = mapVBVD_JB(raw_data);
twix_img = twix{2};
PMU = twix_img.PMUdata;
sum(sum(PMU.EXT))
sum(sum(PMU.raw.EXT.data))

if save_trigger
    disp(['Saving twix into', twix_path]);
    save(twix_path, 'twix');
end
%% April 17 exploring the EVNT
% twix_image2 = twix{1,2};
% pmuEVNT = twix_image2.PMUdata.EVNT;
% pmuTimestamp = twix_image2.image.timestamp;
% 
% stop_time_ms = (pmuTimestamp(end)-pmuTimestamp(1)) * 2.5;
% time_ms = linspace(0, stop_time_ms, length(pmuTimestamp));
% evnt = double(pmuEVNT(1,:));
% 
% pmu_evnt_table = array2table([pmuTimestamp(:)'; time_ms(:)' ; evnt(:)'; ]', 'VariableNames', ...
%     {'PMU_timestamp', 'PMU_time_ms', 'EVNT'});
% pmu_mark_table = pmu_evnt_table;
% pmu_mark_table((pmu_evnt_table.EVNT == 0),:)=[];
%%
twix_image2 = twix{1,2};
pmuEXT = twix_image2.PMUdata.EXT;
pmuTimestamp = twix_image2.image.timestamp;

stop_time_ms = (pmuTimestamp(end)-pmuTimestamp(1)) * 2.5;
time_ms = linspace(0, stop_time_ms, length(pmuTimestamp));
ext1 = double(pmuEXT(1,:));
ext2 = double(pmuEXT(2,:));%remember to convert the data type, otherwise something wrong with rawTimestamp.
%
pmu_ext_table = array2table([pmuTimestamp(:)'; time_ms(:)' ; ext1(:)'; ext2(:)']', 'VariableNames', ...
    {'PMU_timestamp', 'PMU_time_ms', 'EXT1', 'EXT2'});
pmu_mark_table = pmu_ext_table;
pmu_mark_table((pmu_ext_table.EXT1 == 0) & (pmu_ext_table.EXT2 == 0),:)=[];
%%
figure;set(gcf, 'Color', 'w');  % Set figure background to white
nPoints = 4000;
plot(time_ms(1:nPoints), pmu_ext_table.EXT1(1:nPoints), 'r')
hold on;
plot(time_ms(1:nPoints), pmu_ext_table.EXT2(1:nPoints), 'b--')
grid on;
xlabel('time (ms)')
ylabel('pulse amplitude')
legend('EXT1', 'EXT2');
xlim([3300 4000])
%%
% Example data
t_ms = time_ms(:); % Time vector (in seconds)
t_tp = pmuTimestamp(:);
y = pmu_ext_table.EXT2(:); % Example pulse signal with noise

% Find peaks
disp('Make sure Signal Processing Add-on is installed to use findpeaks')
[peaks_amp, locs_time_ms] = findpeaks(y, t_ms); % locs will contain timestamps
[peaks_amp, locs_tp] = findpeaks(y, t_tp); % locs will contain timestamps

% Display results
% disp('Detected Peaks and Corresponding Timestamps:')
% disp(array2table([peaks_amp, locs_time_ms, locs_tp]))
% 
%% Plot
figure; set(gcf, 'Color', 'w');  % Set figure background to white
plot(t_ms, y, 'b'); hold on;
plot(locs_time_ms(:), peaks_amp(:), 'ro', 'MarkerFaceColor', 'r'); % Mark peaks
xlabel('Time (ms)');
ylabel('Amplitude');
xlim([3000 9000])
title('Detected Peaks in Pulse Sequence',  'FontSize', 16);
grid on;
%%
disp('check intervals between the triggers:')
% diff(locs_time_ms)
disp('The interval between the start of MRI and the first trigger:')
locs_time_ms(1)
%-----------------------------------------------------------------

%%
% x = [0 1 2 2 2 1 0];  % Plateau peak at 2 (indices 3-5)
% [pks, locs, widths, proms] = findpeaks(x);
% 
% % Find start and end of plateau (manual method)
% plateau_val = pks(1);  % assuming one peak
% idx = find(x == plateau_val);  % all indices where value is the peak
% 
% % Find connected indices (plateau)
% diffs = diff(idx);
% split_points = [0 find(diffs > 1) length(idx)];
% for i = 1:length(split_points)-1
%     seg = idx(split_points(i)+1 : split_points(i+1));
%     midpoint_idx = seg(floor(end/2)+1);
%     disp(['Midpoint of plateau: ', num2str(midpoint_idx)]);
% end