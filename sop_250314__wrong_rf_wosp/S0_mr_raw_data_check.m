clc; clear;
%% ====================================================
% Author: Yiwei Jia
% Date: March 14, 2025.
% Happy pi day!!
% derived from /Users/cag/Documents/forclone/mapVBVD/twix_process_yj
% --------------------------------------------------------------------------------------------------
% We need to intially inspect MR raw data for
% (1) QC of raw data once a new batch of acquisition is done
% (2) Check the duration of raw data to determine the length of ET mask
% (3) Determine the start of trigger from physio: Now scanner -> trigger
% ====================================================

subject_num = 1;

if subject_num == 1
    %a: rf1: 2 rf pulses
    meas_name_suffix = '_MID00143_FID292458_a1_rf1';
elseif subject_num == 2
    %b: rf2: single rf pulse with 2 parts
    meas_name_suffix = '_MID00151_FID292466_b1_rf2';
else
    %c: rect rf pulses as reference
    meas_name_suffix = '_MID00154_FID292469_c1_rect_rf';
end

meas_name = ['meas', meas_name_suffix];
twix_name = ['twix', meas_name_suffix];
twix_path =  ['/home/debi/yiwei/mreye_dataset/250127_acquisition/', twix_name,'.mat'];
datasetDir = '/home/debi/yiwei/mreye_dataset/250314_libre_seq_rfsp/data/';
raw_data = [datasetDir, meas_name,'.dat'];
save_trigger = false;
%% (1) (2)
inspect_raw_data = true;
if inspect_raw_data
    bmTwix_info(raw_data);
    check_data_timestamp(datasetDir,false);
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
figure;
plot(time_ms(1:1*1000), pmu_ext_table.EXT1(1:1*1000))
%plot(time_ms(1:10*1000), pmu_ext_table.EXT2(1:10*1000))
xlabel('time (ms)')
ylabel('pulse amplitude')
%%
% Example data
t_ms = time_ms(:); % Time vector (in seconds)
t_tp = pmuTimestamp(:);
y = pmu_ext_table.EXT1(:); % Example pulse signal with noise

% Find peaks
[peaks_amp, locs_time_ms] = findpeaks(y, t_ms); % locs will contain timestamps
[peaks_amp, locs_tp] = findpeaks(y, t_tp); % locs will contain timestamps

% Display results
disp('Detected Peaks and Corresponding Timestamps:')
disp(array2table([peaks_amp, locs_time_ms, locs_tp]))

%% Plot
figure;
plot(t_ms, y, 'b'); hold on;
plot(locs_time_ms, peaks_amp, 'ro', 'MarkerFaceColor', 'r'); % Mark peaks
xlabel('Time (ms)');
ylabel('Amplitude');
title('Detected Peaks in Pulse Sequence');
grid on;
%%
disp('check intervals between the triggers:')
% diff(locs_time_ms)
disp('The interval between the start of MRI and the first trigger:')
locs_time_ms(1)
%-----------------------------------------------------------------
% Detected Peaks and Corresponding Timestamps:
%      Var1        Var2          Var3   
%     ______    __________    __________
% 
%     14.875          2424    2.2003e+07
%          1          4928    2.2004e+07
%          1          7448    2.2005e+07
%          1          9952    2.2006e+07
%          1         12448    2.2007e+07
%          1         14952    2.2008e+07
%          1         17448    2.2009e+07
%     14.875         19944     2.201e+07
%         15         22448    2.2011e+07
%          1         24952    2.2012e+07
%          1         27472    2.2013e+07
%          1         29968    2.2014e+07
%          1         32472    2.2015e+07
%          1         34968    2.2016e+07
%-----------------------------------------------------------------
