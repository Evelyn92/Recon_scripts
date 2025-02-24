
%% Set up
clear;clc;



%% Get Path to all to modify Subjects' files
% Prompt user to select Siemens Raw Data file
[fileName, filePath] = uigetfile({'*.dat', '(*.dat) Siemens Rawdatafile'}, ...
    'Select a Subject Rawdatafile', 'MultiSelect', 'off',fileparts(config.foldersPath));
newdatapath = fileparts(fileparts(fileparts(filePath)));

if SINGLEFILEFLAG
    fileList = dir(fullfile(filePath, fileName));
else
    % 27 as BEAT_LIBREoff_BOLD_Cerv.dat is shorter with 27 characters
    fileList=dir(fullfile(newdatapath,'**','*',fileName(end-27:end))); 
    %fileList=dir(fullfile(newdatapath,'**','*BEAT_LIBREoff_BOLD_Cerv.dat'));
    %fileList=dir(fullfile(newdatapath,'**','*BEAT_LIBREoff_BOLD_audio_bis.dat'));
end

% You may also want to output some other variables or just save things
% intermediately as described throughout

% Add paths for subfunctions
addpath(genpath('Dependencies'))

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!!!!!!!!!!! Display the session pid and make sure you log it in the calendar!!!!
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
disp('Please log your pid in the calendar')
feature getpid

%% Start calculations
for iFile=1:length(fileList)
    
    if ~exist(fullfile(fileList(iFile).folder(1:end-5),'5_CorrectionMatrix','CorrectionMatrix.mat'),'file')||OVERWRITEFLAG
        
        
        
        disp(['File: ',num2str(iFile),' of: ',num2str(length(fileList))])
        disp([fileList(iFile).folder,filesep,fileList(iFile).name])
        filepathRawData=[fileList(iFile).folder,filesep,fileList(iFile).name];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Part 1: Read in raw data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        rotTraslFile = dir([fileList(iFile).folder(1:end-5),filesep,'3_SPM_Analysis/2_Vols/*.txt']);
        
        % Take latest if multiple files exist
        if numel(rotTraslFile) > 1
            [~, indexes] = sort([rotTraslFile.datenum]);
            rotTraslFile = rotTraslFile(indexes(end));
        end

        if exist([rotTraslFile.folder,filesep,rotTraslFile.name],'file')
            
            % Load in the precomputed motion estimates
            % load([fileList(iFile).folder(1:end-6),filesep,'3_SPM_Analysis/*.txt'],'rpON1')
            rpON1 = importdata([rotTraslFile.folder,filesep,rotTraslFile.name]);
            % Convert tables in to vectors containing the translational and
            % rotational components
            Translations=[rpON1(:,1)';rpON1(:,2)';rpON1(:,3)'];
            Angles=[rpON1(:,4)';rpON1(:,5)';rpON1(:,6)'];
            clear rpON1

            % Load binning masks so we know which motion components go with which
            % readouts
            % This is the first reconstruction folder: it contains two
            % files imCS and parameters
            load([fileList(iFile).folder,filesep,'reconParams',filesep,'parameters'])
            
            MASK_ONOFF=param.eyeONOFFBinning;
            nONOFF = size(param.eyeONOFFBinning,2);
            nTrials = size(Translations,2)/nONOFF;
            MASK_TRIAL=param.eyeTRIALSBinning;
            clear param
            
            % Reshape motion vectors so we can easily assign the trial or temporal
            % component to the correct k-space readouts
            Translations=reshape(Translations,3,nONOFF,[]);
            Angles=reshape(Angles,3,nONOFF,[]);
            
            % Just permuting some dimensions of the masks to make sure things are
            % consistent. Some of the masks were columns and some were row vectors.
            for i=1:length(MASK_ONOFF)
                MASK_ONOFF{1,i}=MASK_ONOFF{1,i}(:);
            end
            for i=1:length(MASK_TRIAL)
                MASK_TRIAL{i,1}=MASK_TRIAL{i,1}(:);
            end
            
            
            %    Assign the motion components to new vectors that span the entire
            %    acquisition according to the masks
            % We will then use nearest neighbor interpolation to ensure we have a
            % component for every readout
            TRANSLATIONS=zeros(length(MASK_ONOFF{1,1}),3);
            ANGLES=zeros(length(MASK_ONOFF{1,1}),3);
            SAMPLEMASK=false(length(MASK_ONOFF{1,1}),1);
            for iONOFF=1:nONOFF
                for iTRIAL=1:nTrials
                    TRANSLATIONS(MASK_ONOFF{1,iONOFF}&MASK_TRIAL{iTRIAL,1},:)=repmat(Translations(:,iONOFF,iTRIAL)',[sum(MASK_ONOFF{1,iONOFF}&MASK_TRIAL{iTRIAL,1}),1]);
                    SAMPLEMASK(MASK_ONOFF{1,iONOFF}&MASK_TRIAL{iTRIAL,1},1)=true;
                    ANGLES(MASK_ONOFF{1,iONOFF}&MASK_TRIAL{iTRIAL,1},:)=repmat(Angles(:,iONOFF,iTRIAL)',[sum(MASK_ONOFF{1,iONOFF}&MASK_TRIAL{iTRIAL,1}),1]);
                end
            end
            t=linspace(0,1,length(MASK_ONOFF{1,1}));
            TRANSLATIONS=-interp1(t(SAMPLEMASK),TRANSLATIONS(SAMPLEMASK,:),t,'nearest','extrap');
            ANGLES=-interp1(t(SAMPLEMASK),ANGLES(SAMPLEMASK,:),t,'nearest','extrap');
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Begin reading in raw data and header information (twix_obj)
            twix_obj = mapVBVD_JH(filepathRawData);
            if iscell(twix_obj);twix_obj = twix_obj{end};end
            twix_obj.image.flagIgnoreSeg = true; %(essential to read the data correctly)
            
            
            % Store some of the header information related to the sequence
            % parameters
            param.Np     = double(twix_obj.image.NCol);     % Number of readout point per spoke
            if twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{13}==2
                param.Nshot  = floor(double(twix_obj.image.NSeg)/twix_obj.hdr.MeasYaps.sWipMemBlock.alFree{14});     % Number of heartbeats
            else
                param.Nshot  = floor(double(twix_obj.image.NSeg));
            end
            param.Nseg   = double(twix_obj.image.NLin/param.Nshot); % Number of segments per heartbeat
            param.Nx        = param.Np/2;
            param.Ny        = param.Np/2;
            param.Nz        = param.Np/2;
            
            
            
            
            
            [kx, ky, kz] = computePhyllotaxis (param.Np, param.Nseg, param.Nshot, true, false, false);
            % For every readout, rotate the phylotaxis trajectory by the
            % corresponding rorational component of the motion estimations
            TRAJ=zeros(3,param.Np,param.Nseg*param.Nshot);
            for iSpoke=1:size(ANGLES,1)
                TRAJ(:,:,iSpoke) = my3DROT([col(kx(:,iSpoke)),col(ky(:,iSpoke)),col(kz(:,iSpoke))]',ANGLES(iSpoke,:),true);
            end
            newkx=reshape(TRAJ(1,:,:),[param.Np,param.Nseg,param.Nshot]);
            newky=reshape(TRAJ(2,:,:),[param.Np,param.Nseg,param.Nshot]);
            newkz=reshape(TRAJ(3,:,:),[param.Np,param.Nseg,param.Nshot]);
            
            PHASERAMP=single(exp(2*pi*1i*(kx(:,:).*TRANSLATIONS(:,1)'+ky(:,:).*TRANSLATIONS(:,2)'+kz(:,:).*TRANSLATIONS(:,3)')));
            
            
            if ~exist([fileList(iFile).folder(1:end-5),filesep,'5_CorrectionMatrix',filesep],'dir')
                mkdir([fileList(iFile).folder(1:end-5),filesep,'5_CorrectionMatrix',filesep])
            end

            % Get the current date and time
            currentDateTime = datetime('now', 'Format', 'yyyyMMdd_HHmm');
    
            % Create a formatted string with the desired format
            formattedString = sprintf('%04d-%02d-%02d_%02d-%02d', ...
            year(currentDateTime), month(currentDateTime), day(currentDateTime), ...
            hour(currentDateTime), minute(currentDateTime));

            % Create filename
            mocoFileName = ['CorrectionMatrix_',formattedString];
            
            save([fileList(iFile).folder(1:end-5),filesep,'5_CorrectionMatrix',filesep,mocoFileName],'PHASERAMP','newkx','newky','newkz','-v7.3')
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Optional quick static (using all lines except SI) recon to make sure MOCO worked
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if RECONFLAG
                
                % Read in the actual KSpace data
                KSPACE = twix_obj.image.unsorted();KSPACE = permute( KSPACE, [1 3 4 2] ); % Dimensions of rawData are now {'Col'  'Lin'  'Set'  'Cha'}
                
                % Original trajectory
                [kx, ky, kz] = computePhyllotaxis (param.Np, param.Nseg, param.Nshot, true, false, false);
                Kstatic   = cat(4, kx, ky, kz);
                Kstatic   = single(reshape(Kstatic,[param.Np, param.Nseg*param.Nshot, 1, 3]));
                
                % Corrected trajectory trajectory
                newKstatic   = cat(4, newkx, newky, newkz);
                newKstatic   = single(reshape(newKstatic,[param.Np, param.Nseg*param.Nshot, 1, 3]));
                
                % Create a mask that will remove the SI projections from the data
                % acquired at the beginning of each shot as well as the first 20 shots
                % before steady stade is reached
                mask=true(param.Nseg*param.Nshot,1);
                mask(1:param.Nseg:end)=false;
                mask(1:param.Nseg*20)=false;
                
                % Calculate the density compensation weights for the non cartesian
                % reconstruction
                Wstatic = densityCompensationForUniform3DRadial(param.Np, param.Nseg*param.Nshot, true);
                Wstatic = repmat(Wstatic(:), [1, param.Nseg*param.Nshot]);
                
                % Uncorrected reconstruction
                NUFFT = gpuNUFFT((reshape(Kstatic(:,mask ,1,:),[param.Np*sum(mask),3]))',(col(Wstatic(:,mask))).^2,1.5,3,8,[param.Nx,param.Ny,param.Nz],[],true);%clear Kstatic clear Wstatic
                I0=sos(NUFFT'*(reshape(KSPACE(:,mask ,1,:),[param.Np*sum(mask),size(KSPACE,4)])));
                
                % Corrected reconstruction
                NUFFT = gpuNUFFT((reshape(newKstatic(:,mask ,1,:),[param.Np*sum(mask),3]))',(col(Wstatic(:,mask))).^2,1.5,3,8,[param.Nx,param.Ny,param.Nz],[],true);%clear Kstatic clear Wstatic
                IMOCO=sos(NUFFT'*(reshape(PHASERAMP(:,mask).*KSPACE(:,mask ,1,:),[param.Np*sum(mask),size(KSPACE,4)])));
                
                % Uncomment if you want to compare the reconstructions
                    im5D(cat(4,I0,IMOCO))
            end
            
        end
    end
end