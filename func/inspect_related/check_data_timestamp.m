function check_data_timestamp(databasedir,inspect)


[rawDataName, rawDataDir, ~] = uigetfile( ...
    { '*.dat','Siemens raw data file (*.dat)'}, ...
       'Pick a file', ...
       'MultiSelect', 'off',databasedir);

    if rawDataName == 0
        warning('No file selected');
        return;
    end

    filepathRawData = fullfile(rawDataDir, rawDataName);
%--------------------------------------------------------------------------
% Read the PMUTimeStamp triggered from external trigger
%--------------------------------------------------------------------------  
    costTime = 2.5;
    param.batchParam.rawDataName    = rawDataName;
    param.batchParam.rawDataDir     = rawDataDir;

    [ twix_obj, param ] = dataSelectionAndLoading( databasedir, param );
    % Load the PMUTime and TimeStamp, shift the TimeStamp to the beginning
    % of 0
    PMUTimeStamp    = double( twix_obj.image.pmutime );
    TimeStamp       = double( twix_obj.image.timestamp );

    %
    TimeStamp       = TimeStamp - min(TimeStamp);
    % Do not forget to scale the times by costTime (Setting from Siemens)
    PMUTimeStamp_ms = PMUTimeStamp * costTime;
    PMUTimeStamp_s  = PMUTimeStamp_ms / 1000;
    TimeStamp_ms    = TimeStamp * costTime;
    TimeStamp_s     = TimeStamp_ms / 1000;
    % Save all the param to struct "param"
    param.PMUTimeStamp_ms   = PMUTimeStamp_ms;
    param.PMUTimeStamp_s    = PMUTimeStamp_s;
    param.TimeStamp_ms      = TimeStamp_ms;
    param.TimeStamp_s       = TimeStamp_s;


if inspect
    figure;
    plot(param.TimeStamp_ms, param.PMUTimeStamp_ms);
    % xlim([0 1e4]);  % Set x-axis limits
    % ylim([0 180]); % Set y-axis limits
    xlabel('TimeStamp (ms)');
    ylabel('PMU TimeStamp (ms)');
    title('Time vs. PMU Stamp');
end


Timediff = TimeStamp_ms(end) - TimeStamp_ms(1);

disp('In protocol, we set the T1 LIBRE duration: 362(sec)*1000(freq)')
disp('And we set the T1 VIBE duration: 357(sec)*1000(freq)')
disp(['However, here the duration of the rawdata is: ', num2str(Timediff), ...
    ' ms with data points:', num2str(length(TimeStamp_s)) ]);

end
