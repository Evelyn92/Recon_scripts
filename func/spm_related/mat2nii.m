% Define the path to your .mat file and the output NIfTI file
x_woBinning = load('/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub002/T1_LIBRE_woBinning/output/x_nIter3.mat');
x_woBinning = x_woBinning.x;
volume_data = x_woBinning{1};

% Replace with your desired output NIfTI file path
nifti_file = '/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub002/T1_LIBRE_woBinning/output/x_nIter3.nii';       




% Ensure the volume_data is in the correct format (3D or 4D matrix)
if ndims(volume_data) < 3
    error('The data in the .mat file must be at least 3D.');
end

%%
% Define NIfTI metadata (optional but recommended for completeness)
% You can adjust these properties according to your needs.
nii_hdr = struct;  % Create default NIfTI header
nii_hdr.ImageSize = size(volume_data);
nii_hdr.PixelDimensions = [0.5 0.5 0.5];  % Adjust these values if needed

% Write the NIfTI file
% niftiwrite(volume_data, nifti_file, nii_hdr);
niftiwrite(volume_data, nifti_file);
disp(['Data has been saved as a NIfTI file: ', nifti_file]);
%%
sub001_woBinning = '/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub001/T1_LIBRE_woBinning/output/x_nIter3.nii';
sub002_woBinning ='/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub002/T1_LIBRE_woBinning/output/x_nIter3.nii';

sub001_woBinning_png = '/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub001/T1_LIBRE_woBinning/output/x_nIter3_png_nii';
sub002_woBinning_png = '/media/sinf/1,0 TB Disk/240825_recon_for_poster/Sub002/T1_LIBRE_woBinning/output/x_nIter3_png_nii';

% chop_volume(sub001_woBinning, sub001_woBinning_png)
chop_volume(sub002_woBinning, sub002_woBinning_png)