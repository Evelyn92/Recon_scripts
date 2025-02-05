function co_register(source_image,ref_image)
%CO_REGISTER 
%  co_register the image source to the image ref
%  There is no output since it is an in-place operation
disp(['The path to source image: ', source_image])
disp(['The path to reference image: ', ref_image])

matlabbatch{1}.spm.spatial.coreg.estimate.ref = {[ref_image,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.source = {[source_image,',1']};
matlabbatch{1}.spm.spatial.coreg.estimate.other = {''}; % No additional images
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi'; % Normalized Mutual Information for alignment
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.sep = [4 2]; % Sampling at 4mm and 2mm
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001]; % Tolerance for the alignment
matlabbatch{1}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7]; % Smoothing kernel for coregistration

% Run the batch
spm_jobman('run', matlabbatch);


end

