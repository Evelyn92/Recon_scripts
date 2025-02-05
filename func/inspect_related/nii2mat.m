function nii2mat(base_folder, equal)
    
    [niiName, niiDataDir, ~] = uigetfile( ...
        { '*.nii','nii file (*.nii)'}, ...
           'Pick a nii file to be converted', ...
           'MultiSelect', 'off', base_folder);

    if niiName == 0
        warning('No nii file selected');
        return;
    end

    filepathniiData = fullfile(niiDataDir, niiName);
    disp(filepathniiData);
    [~, niiName_only, ~] = fileparts(niiName);

 
    % Load the NIfTI file using SPM
    nii = spm_vol(filepathniiData);  % Load NIfTI file (can be .nii or .nii.gz)
    data = spm_read_vols(nii);       % Read the volumetric data
    data = rot_norm_equal(data, equal);
    output_folder = [niiDataDir, '/', 'vibe_mat/'];
    if ~exist(output_folder, 'dir')
        mkdir(output_folder);
    end

    % Save the data to a .mat file
    matSavePath = [output_folder,'/', niiName_only, '.mat'];
    save(matSavePath, 'data');
    disp(['The converted mat file is saved here: ', matSavePath]);

end