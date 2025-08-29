
%% Unzip Eye atlas .nii.gz files to .nii

gunzip('/Users/cag/Documents/Dataset/A-Eye-altas/eye_atlas/female/*.nii.gz');
gunzip('/Users/cag/Documents/Dataset/A-Eye-altas/eye_atlas/male/*.nii.gz');

%% Convert MAT to Nii format

%% Add SPM to MATLAB path and start SPM
addpath('/Users/cag/Documents/forclone/spm'); spm('Defaults','fMRI'); spm_jobman('initcfg');

% Define paths
atlas_template = '/Users/cag/Documents/Dataset/A-Eye-altas/eye_atlas/male/template.nii';     % E.g. MR‑Eye template
maxprob = '/Users/cag/Documents/Dataset/A-Eye-altas/eye_atlas/male/max_prob_map.nii';         % Max‑prob map with labeled ROIs
roi_labels = 1:9;  % Assuming 9 eye structures in A‑Eye atlas

% List your subjects
subjects = {'subj01','subj02'};  % replace with your IDs
data_dir = '/path/to/your/data';

% for i = 1:length(subjects)
for i = 1
    subj = subjects{i};
    subj_t1 = fullfile(data_dir, subj, 'anat', [subj '_T1w.nii']);

    % Step 1: Normalize (Estimate & Write) to atlas template space
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {subj_t1};
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = {subj_t1};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.template = {atlas_template};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smosrc = 8;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.smoref = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.regtype = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.cutoff = 25;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.nits = 16;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = 1;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.preserve = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.bb = [-90 -126 -72; 90 90 108];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.vox = [1 1 1];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.normalise.estwrite.roptions.prefix = 'w';

    spm_jobman('run', matlabbatch);
    clear matlabbatch;

    % Step 2: Load normalized subject and extract signal per ROI
    wsubj_t1 = fullfile(data_dir, subj, 'anat', ['w' subj '_T1w.nii']);
    subj_vol = spm_vol(wsubj_t1);
    subj_img = spm_read_vols(subj_vol);
    
    roi_vol = spm_vol(maxprob);
    roi_img = spm_read_vols(roi_vol);

    noise_mask = (roi_img == 0);  
    sigma_noise = std(subj_img(noise_mask));

    for lbl = roi_labels
        mask = (roi_img == lbl);
        roi_signal = subj_img(mask);
        mean_sig = mean(roi_signal);
        vol_size = sum(mask(:));
        results(i,lbl).mean_signal = mean_sig;
        results(i,lbl).volume = vol_size;
        results(i,lbl).SNR = mean_sig / sigma_noise;
    end
end

% Step 3: Compute CNR between ROIs for each subject
for i = 1:length(subjects)
    for a = roi_labels
        for b = roi_labels
            if a < b
                results(i,a).CNRto(b) = abs(results(i,a).mean_signal - results(i,b).mean_signal) / sigma_noise;
            end
        end
    end
end

% Display summary
disp(results);