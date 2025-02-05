function reslicing(source_image,ref_image)
%% Step 1: Coregister and Reslice Source Image B to match Reference Image A

matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[ref_image,',1']}; % Reference Image A
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[source_image,',1']}; % Source Image B
matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''}; % No other images to coregister
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi'; % Normalized mutual information for cost function
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2]; % Sampling at 4mm and 2mm steps
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; % Tolerances
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7]; % Smoothing kernel for coregistration
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4; % 4th-degree spline interpolation
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0]; % No wrapping in x, y, z directions
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0; % No masking
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'resliced_'; % Prefix for resliced image

% Run the batch
spm_jobman('run', matlabbatch);
disp('The path to the resliced image:');
% Extract directory and filename
[dir_path, file_name, ext] = fileparts(source_image);
% Generate new filename with prefix
new_file_name = ['resliced_', file_name, ext];
% Construct new path
resliced_path = fullfile(dir_path, new_file_name);
% Display the new path
disp(resliced_path);

end