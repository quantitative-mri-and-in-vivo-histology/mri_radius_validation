function compute_effective_radius_power_law_ratio_dwi( ...
    powder_average_file, ...
    bval_file, ...
    output_file, ...
    delta, ...
    Delta, ...
    D0, ...
    options)
% Estimates effective radius using a power-law ratio of powder-averaged
% signals over two diffusion shells. The approach was described by 
% Pizzolato et. al (https://doi.org/10.1016/j.media.2023.102767)
%
% INPUTS:
%   input_file  - (string) Path to the input 4D NIfTI file containing 
%                 powder-averaged signals.
%   bval_file   - (string) Path to the corresponding b-value file.
%   output_file - (string) Path where the computed effective radius map 
%                 should be saved.
%   delta       - (double) Diffusion gradient pulse duration (in ms).
%   Delta       - (double) Diffusion gradient separation time (in ms).
%   D0          - (double) Free diffusion coefficient (in micron²/ms).
%
% OPTIONS (Name-Value Pairs):
%   options.mask_file    - (string, default = "") Path to an optional 
%                          binary mask file. If provided, computations 
%                          are performed only within the masked region.
%   options.min_bval     - (numeric, default = 0) Minimum b-value to 
%                          include in calculations (in ms/micron²).
%   options.max_bval     - (numeric, default = inf) Maximum b-value to 
%                          include in calculations (in ms/micron²).
%
% OUTPUTS:
%   The function writes the computed effective radius map as a NIfTI file 
%   to 'output_file'. If 'options.d_perp_file' is specified, the computed 
%   D_perp map is also saved as a NIfTI file.
%
% USAGE:
%       compute_effective_radius_power_law_ratio_dwi("dwi.nii.gz", ... 
%       'dwi.bval', 'r_eff.nii.gz', 15, 29.25, 2.07, ...
%       'mask_file', 'mask.nii.gz');

    arguments
        powder_average_file string
        bval_file string
        output_file string
        delta double
        Delta double
        D0 double   
        options.mask_file (1,1) string = "" 
        options.min_bval (1,1) {mustBeNumeric} = 0
        options.max_bval (1,1) {mustBeNumeric} = inf
    end

    bvals = importdata(bval_file); 
    bval_idx = (bvals >= options.min_bval) & (bvals <= options.max_bval);
    bvals = bvals(bval_idx)/1000;

    powder_average_image = niftiread(powder_average_file);
    powder_average_image = double(powder_average_image);
    powder_average_image = powder_average_image(:,:,:,bval_idx);
    invalid_idx = powder_average_image < 0 ...
        | isinf(powder_average_image) ...
        | isnan(powder_average_image);
    powder_average_image(invalid_idx) = 0;

    if options.mask_file ~= ""
        mask = niftiread(options.mask_file);
        mask_rep = repmat(mask, 1, 1, 1, size(powder_average_image,4));
        powder_average_image(~mask_rep) = 0;
    end
    
    invalid_idx = zeros(size(powder_average_image,1:3));
    for i = 1:size(powder_average_image,4)
        invalid_idx = invalid_idx ...
            | (powder_average_image(:,:,:,i) == 0);
    end

    [b_high, b_high_idx] = max(bvals);
    [b_low, b_low_idx] = min(bvals);
    powder_average_high_shell = powder_average_image(:,:,:,b_high_idx);
    powder_average_low_shell = powder_average_image(:,:,:,b_low_idx);

    D_perp = log(powder_average_low_shell./powder_average_high_shell ...
        .*sqrt(b_low./b_high))./(b_high-b_low);
    r_eff_map = (48/7.*delta.*(Delta-delta/3).*D0.*D_perp).^(1/4);
    r_eff_map(invalid_idx | D_perp < 0) = nan;

    use_compression = false;
    if endsWith(output_file, ".nii.gz")
        use_compression = true;
        output_file = erase(output_file, ".nii.gz");
    elseif endsWith(output_file, ".nii")
        output_file = erase(output_file, ".nii");
    end
    nifti_info = niftiinfo(powder_average_file);
    nifti_info.ImageSize = nifti_info.ImageSize(1:ndims(r_eff_map));
    nifti_info.PixelDimensions = ...
        nifti_info.PixelDimensions(1:ndims(r_eff_map));
    nifti_info.Datatype = class(r_eff_map);
    niftiwrite(r_eff_map, ...
        output_file, ...
        nifti_info, ...
        'Compressed', ...
        use_compression);
end
