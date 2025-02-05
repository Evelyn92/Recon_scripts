function [data_rot] = rot_norm_equal(data, equal)
%ROT_NORM_EQUAL Summary of this function goes here
%   Detailed explanation goes here
data_rot = rand(size(data));

% Rotate the (x, y) plane by 90 degrees for each slice along z
for k = 1:size(data, 3)
    slice = data(:, :, k);
    slice = rot90(slice, -1);
    slice = abs(slice);
    slice = mat2gray(slice);
  
    if equal
        min_val = 0.01;
        max_val = 0.5;
        slice = equal_func(slice, min_val, max_val);
    end
    
    data_rot(:, :, k) = slice;
end


end

