%--------------------------------------------------------------------------
% custom change to import mapVBVD: not sure its correct, no time to test it.
% Just need to be sure you are able to acces the mapVBVD_JH function for example.
%--------------------------------------------------------------------------
addpath('C:\yiwei\Recon_fork\MR_EYE_RawData_Binning_fork\code\MATLAB\ReadRawDataSiemens\mapVBVD')
%%
%--------------------------------------------------------------------------
% Original code from the function BinningEYE.m
%--------------------------------------------------------------------------
    param.basedir = 'C:\yiwei\1_Pilot_MREye_Data\Sub001\230928_anatomical_MREYE_study\MR_EYE_Subj01\RawData';
    basedir = param.basedir;
    param.batchParam = [];
%%
% Prompt user to select Siemens Raw Data file
%--------------------------------------------------------------------------

    [rawDataName, rawDataDir, ~] = uigetfile( ...
    { '*.dat','Siemens raw data file (*.dat)'}, ...
       'Pick a file', ...
       'MultiSelect', 'off',basedir);

    if rawDataName == 0
        warning('No file selected');
        return;
    end

    filepathRawData = fullfile(rawDataDir, rawDataName);
%%
%--------------------------------------------------------------------------
% Initialize saving directory
%--------------------------------------------------------------------------    
  % Create saving directory
    [~,name,~] = fileparts(rawDataName);
%%
  % Add number of bin to saving directory name
    name = sprintf('BinningEYE_%s',name);
    
  % Define directory name to save the results of the ICA analysis
    param.savedir = fullfile(rawDataDir,name);
    
  % Create directory
    if exist(param.savedir,'dir') ~= 7
        mkdir(param.savedir);
    end
%%
%--------------------------------------------------------------------------
% Prompt user to enter the number of desired bins
%-------------------------------------------------------------------------- 
% How to select the number of bins?
    prompt        = {'Enter the number of bins'};
    name          = '#Bins';
    numlines      = 1;
    defaultanswer = {'4'};
    answer        = inputdlg(prompt,name,numlines,...
                           defaultanswer);
    
  % Check if the user selected cancel
    if isempty(answer)
        warning('The user selected cancel');
        return;
    end
    
  % Convert string answer to number
    nbins = str2double(answer{1});
    param.nBins = nbins;
    
%%
%--------------------------------------------------------------------------
% Read the PMUTimeStamp triggered from external trigger
%--------------------------------------------------------------------------  
    costTime = 2.5;
    param.batchParam.rawDataName    = rawDataName;
    param.batchParam.rawDataDir     = rawDataDir;

    [ twix_obj, param ] = dataSelectionAndLoading( basedir, param );
    % Load the PMUTime and TimeStamp, shift the TimeStamp to the beginning
    % of 0
    PMUTimeStamp    = double( twix_obj.image.pmutime );
    TimeStamp       = double( twix_obj.image.timestamp );
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
%% Inspect the visualization of PMUTime and TimeStamp
    figure;
    plot(param.TimeStamp_ms, param.PMUTimeStamp_ms);
    % xlim([0 1e4]);  % Set x-axis limits
    % ylim([0 180]); % Set y-axis limits
    xlabel('TimeStamp (ms)');
    ylabel('PMU TimeStamp (ms)');
    title('Time vs. PMU Stamp');
%%
    [ triggerPeaks_ms, triggerIntervals_ms ] = findpeaks_PMUTimeStamp( param.PMUTimeStamp_ms, param.TimeStamp_ms );
    % length(triggerPeaks_ms)=106 length(triggerIntervals_ms)=104
    triggerPeaks_ms = [ -1 triggerPeaks_ms ];
    triggerPeaks_s      = triggerPeaks_ms       / 1000;
%     triggerIntervals_s  = triggerIntervals_ms   / 1000;
    nTrig = length( triggerPeaks_ms );
% %%   Explore the syncboxMatrix
% %--------------------------------------------------------------------------
% % Read the sync-box data
% %--------------------------------------------------------------------------
% 
%     % Prompt user to select respiratory binning directory
%     [ syncboxDataName, syncboxDataDir, ~] = uigetfile( ...
%     { '*.mat','Syncbox data file (*.mat)'}, ...
%        'Pick a file', ...
%        'MultiSelect', 'off',rawDataDir);
% 
%     if syncboxDataName == 0
%         error('No syncbox data file were selected');
%     end
% 
%     syncboxDataPath = fullfile(syncboxDataDir,syncboxDataName);
% 
%     struct = load(syncboxDataPath);
%     syncboxMatrix = struct.data;
%     % what kind of information is saved here at syncboxMatrix?
%     % syncboxMatrix (526640           7)
% %--------------------------------------------------------------------------
% % 
% %--------------------------------------------------------------------------
% 
%     x_mean = syncboxMatrix( :, 5 );
%     y_mean = syncboxMatrix( :, 6 );
% 
%     x_mean = x_mean( x_mean~=0  );
%     y_mean = y_mean( y_mean~=0  );
%     stimuli_pos = [x_mean y_mean];
%     % Size of stimuli_pos is equal to trigger peaks
%     % So the dimension 5 and 6 means the positions of stimuli.
%     % Check the stimuli here
%     scatter(x_mean, y_mean)
%    % 
% 
%     indexVec = syncboxMatrix( :, 3 );
%     indexVec = indexVec( indexVec ~= 0 );
%     %
%     % visualization
%     % Sample data
% 
%     % Generate a color map with 17 unique colors
%     colors = hsv(17); % You can also use other colormap functions like hsv(17), jet(17), etc.
% 
%     % Create a new figure
%     figure;
%     hold on;
% 
%     % Plot each class with a unique color
%     for class = 1:17
%         % Find indices of points belonging to the current class
%         idx = indexVec == class;
% 
%         % Plot the points for the current class
%         scatter(stimuli_pos(idx, 1), stimuli_pos(idx, 2), 50, colors(class, :), 'filled');
%     end
% 
%     % Add labels and title
%     xlabel('Dimension 1');
%     ylabel('Dimension 2');
%     title('Scatter Plot of 106 Points Divided into 17 Classes');
%     legend(arrayfun(@(x) sprintf('Class %d', x), 1:17, 'UniformOutput', false), 'Location', 'BestOutside');
% 
%     % Adjust plot
%     hold off;
%     grid on;
% 
% %     indexVec = indexVec(3:end);
% %%    
%     nPres = length(indexVec);
% 
%     indexVec_x = [];
%     for k = 1:nPres
%         elem = indexVec(k);
%         if elem  ~= 1
%             res = mod((elem-1),4);
%             if res == 0
%                 res = 4;
%             end
%         else
%             % Center
%             res = 0; 
%         end
%         indexVec_x = [ indexVec_x res ];
%     end
% 
%     figure, plot( indexVec_x, '*')
%         ylim([0 5])
%         xlim([0 nPres+1])
%         xlabel('# Presentation')
%         ylabel('Bin')
%         title('Binning along X')
% 
%     indexVec_y = [];
%     for k = 1:nPres
%         elem = indexVec(k);
%         if elem  ~= 1
%             ceilVal = ceil((elem-1)/4);
%         else
%             ceilVal = 0;
%         end
%         indexVec_y = [ indexVec_y ceilVal ];
%     end
% 
%     figure, plot( indexVec_y, 'o')
%         ylim([0 5])
%         title('Binning along Y')
% 
%     figure, plot( indexVec_x, '*')
%     hold on, plot( indexVec_y, 'o')
%         ylim([0 5])
%         xlim([0 nPres+1])
%         xlim([0 nPres+1])
%         xlabel('# Presentation')
%         ylabel('Bin')
%         title('Binning along X and Y')
%         legend('X','Y')
% 
% 
%     meanTR = mean( diff(param.TimeStamp_ms) );
%     figure, plot(diff(param.TimeStamp_ms),'.-')
%         title(['Mean TR = ' num2str(meanTR) ' ms'])
   
%%     
%--------------------------------------------------------------------------
% Assign a bin to each line (X)
%--------------------------------------------------------------------------
    NLin    = length(param.PMUTimeStamp_ms);
    binMask = cell(nbins,1);
    binCnt  = zeros(nbins,1);
    % Here nbins is 4, can I set it to 16 in this case?
    binMaskMatrix = zeros([NLin,nbins]);
    % number of timestamps, number of bins
    trigCnt = 1;
    blockCounter = 0;
    
    for k = 1:NLin
        
        if trigCnt < nTrig
            timeTrig        = triggerPeaks_ms( trigCnt );
            timeTrigPlus1   = triggerPeaks_ms( trigCnt+1 );
            timeElem        = TimeStamp_ms( k );
        end
            
        if timeElem >= timeTrigPlus1
            trigCnt = trigCnt + 1;
            blockCounter = 1;
        else
            % blockCounter ranges from 1 to the upper limit of current
            % interval.
            blockCounter = blockCounter + 1;
        end
        
        % In the case below, they adapt the method only for nbin=4
        % Thus, they generate binning mask according to x axis
        if trigCnt <= nPres
            idx = trigCnt;
            idx_x = indexVec_x( idx );
            % preserve the data when the eye is not gazing at the center
            if idx_x ~= 0
                % why here blockCounter mod 22?
                if mod(k,22) ~=1
                    if blockCounter > 22
                      
                        binMaskMatrix( k, idx_x ) = 1;
                    end
                end
            end
        end
        
    end
      % I guess it is based on the assumption that
                        % the eye becomes stable after 22 timestamp
                        % intervals.
    
    for k = 1:nbins
        binMask{k} = binMaskMatrix(:,k);
    end
    
    param.binMask = binMask;
    
    
%%
%--------------------------------------------------------------------------
% Saving data
%--------------------------------------------------------------------------    
  % Save data
    save(fullfile(param.savedir,'EyeBinning_X_inspect.mat'),'param');
  
    
    
%%   
%--------------------------------------------------------------------------
% Assign a bin to each line (Y)
%--------------------------------------------------------------------------
    NLin    = length(param.PMUTimeStamp_ms);
    binMask = cell(nbins,1);
    binCnt  = zeros(nbins,1);
    
    binMaskMatrix = zeros([NLin,nbins]);
    
    trigCnt = 1;
    blockCounter = 0;
    
    for k = 1:NLin
        
        if trigCnt < nTrig
            timeTrig        = triggerPeaks_ms( trigCnt );
            timeTrigPlus1   = triggerPeaks_ms( trigCnt+1 );
            timeElem        = TimeStamp_ms( k );
        end
            
        if timeElem >= timeTrigPlus1
            trigCnt = trigCnt + 1;
            blockCounter = 1;
        else
            blockCounter = blockCounter + 1;
        end
        
        if trigCnt <= nPres
            idx = trigCnt;
            idx_y = indexVec_y( idx );
            if idx_y ~= 0
                if mod(k,22) ~=1  
                    if blockCounter > 22
                        binMaskMatrix( k, idx_y ) = 1;
                    end
                end
            end
        end
        
    end
    
    for k = 1:nbins
        binMask{k} = binMaskMatrix(:,k);
    end
    
    param.binMask = binMask;
    
    

%--------------------------------------------------------------------------
%% Saving data
%--------------------------------------------------------------------------    
  % Save data
    save(fullfile(param.savedir,'EyeBinning_Y_inspect.mat'),'param');
    
    
    
    