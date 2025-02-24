function compute_powder_average_dwi( ...
    input_file, ...
    bvec_file, ...
    bval_file, ...
    output_file, ...
    options) 
% Computes the powder-average mean of a DWI dataset.
%
%   Computes the spherical mean of diffusion-weighted signals across 
%   different b-value shells and saves the result in a NIfTI file.
%
% inputs:
%   input_file  - (string) Path to the input 4D DWI NIfTI file.
%   bvec_file   - (string) Path to the b-vector text file.
%   bval_file   - (string) Path to the b-value text file.
%   output_file - (string) Path to save the output NIfTI file containing 
%                 the spherical mean DWI.
%
% optional parameters (options structure):
%   'mask_file', mask_file        - (string, default="") Path to a binary 
%                                    mask NIfTI file.
%   'noisemap_file', noisemap_file - (string, default="") Path to a noise 
%                                     map NIfTI file for Rician noise 
%                                     correction.
%   'sph_order', sph_order   - (integer, default=6) Spherical harmonics
%                                   order
%
% usage:
%   % Compute powder-average with a mask and normalization
%   compute_powder_average_dwi('dwi.nii.gz', 'dwi.bvec', 'dwi.bval', ...
%       'dwi_sph_mean.nii.gz', 'mask_file', 'brain_mask.nii.gz', ...
%       'normalize', true);
%

    arguments
        input_file string
        bvec_file string
        bval_file string
        output_file string
        options.mask_file (1,1) string = "" 
        options.noisemap_file (1,1) string = ""
        options.sph_order (1,1) {mustBeInteger, mustBeNumeric} = 6
    end

    input_image = niftiread(input_file); 
    input_image = double(input_image);

    mask = true(size(input_image,1:3));
    if options.mask_file ~= ""
        mask = niftiread(options.mask_file); 
        mask = logical(mask);
    end

    if options.noisemap_file ~= ""
        noisemap = niftiread(options.noisemap_file); 
    end

    mask_idx = find(mask>0);
    bvec = importdata(bvec_file);
    bval = importdata(bval_file);
    bval_unique =  unique(round(bval./10).*10);
    bval_unique_nonzero = bval_unique(bval_unique > 0);
    powder_average_image = zeros([size(mask), length(bval_unique_nonzero)]);

    for i = 1:length(bval_unique_nonzero)
        bval_idx = abs(bval - bval_unique_nonzero(i)) < 500;
        bval_shell_image = input_image(:,:,:,bval_idx);
        
        bval_shell_values = zeros(length(mask_idx), size(bval_shell_image,4));
        for j = 1:length(mask_idx)
            [idx1, idx2, idx3] = ind2sub(size(mask), mask_idx(j));
            bval_shell_values(j,:) = bval_shell_image(idx1,idx2,idx3,:);
        end
        
        pa_options = {};
        if options.noisemap_file ~= ""
            pa_options = {'sigma', noisemap(mask_idx)'};
        end

        powder_average_values = compute_powder_average( ...
            bval_shell_values', ...
            bvec(:,bval_idx)', ...
            options.sph_order, ...
            pa_options{:});

        powder_average_image_temp = zeros(size(mask));
        powder_average_image_temp(mask_idx) = powder_average_values;
        powder_average_image(:,:,:,i) = powder_average_image_temp;
    end
    
    % write powder average to nifti file
    use_compression = false;
    if endsWith(output_file, ".nii.gz")
        use_compression = true;
        output_file = erase(output_file, ".nii.gz");
        bval_file = erase(output_file, ".nii.gz") + ".bval";
    elseif endsWith(image_file, ".nii")
        output_file = erase(output_file, ".nii");
        bval_file = erase(output_file, ".nii") + ".bval";
    end
    nifti_info = niftiinfo(input_file);
    nifti_info.ImageSize = size(powder_average_image);
    nifti_info.PixelDimensions = nifti_info.PixelDimensions(1:ndims(powder_average_image));
    nifti_info.Datatype = class(powder_average_image);
    niftiwrite(powder_average_image, output_file, nifti_info, 'Compressed', use_compression);

    % write powder average bval file
    writematrix(bval_unique_nonzero, bval_file, 'Delimiter', ' ', 'FileType', 'text');
end

