function slice_redistributed = equal_func(slice_norm, min_val, max_val)


slice_norm = mat2gray(slice_norm); % Normalizes the image

% Clip the intensity values
slice_I_clipped = slice_norm;
slice_I_clipped(slice_norm < min_val) = min_val;
slice_I_clipped(slice_norm > max_val) = max_val;

% Normalize to [0, 1] after clipping
slice_redistributed = (slice_I_clipped - min(slice_I_clipped(:))) / (max(slice_I_clipped(:)) - min(slice_I_clipped(:)));


end