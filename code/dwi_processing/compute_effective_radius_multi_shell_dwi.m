function compute_effective_radius_multi_shell_dwi( ...
    powder_average_file, ...
    bval_file, ...
    output_file, ...
    delta, ...
    Delta, ...
    D0, ...
    model, ...
    options)
% Estimates the effective axon radius from multi-shell dMRI data.
%
%   This function computes the effective axon radius (r_eff) voxel-wise using 
%   multi-shell diffusion-weighted MRI (dMRI) data via non-linear fitting
%   to the powder-averaged signal.
%
% INPUTS:
%   powder_average_file - (string) Path to the 4D NIfTI file containing 
%                           powder-average per shell.
%   bval_file           - (string) Path to the file containing b-values 
%                           per powder-average (shell).
%   output_file         - (string) Path to save the estimated effective 
%                           radius map.
%   delta               - (double) Diffusion gradient pulse duration (ms).
%   Delta               - (double) Diffusion gradient separation (ms).
%   D0                  - (double) Free diffusion coefficient 
%                           (micrometer²/ms).
%   model               - (string) Name of the model used for r_eff 
%                           estimation ('Neuman' or 'VanGelderen')
%
% OPTIONAL PARAMETERS (specified in `options` struct):
%   options.mask_file     - (string) Path to a binary mask file (NIfTI).
%                           If provided, computations are performed only 
%                           within the mask.
%   options.beta_file     - (string) Path to save the estimated signal 
%                           scaling factor (beta) 
%   options.min_bval      - (numeric) Minimum b-value threshold to select
%                           powder-average shells (default: 0).
%   options.max_bval      - (numeric) b-value threshold to select
%                           powder-average shells (default: inf).
%
% OUTPUT:
%   A NIfTI file containing voxel-wise estimates of the effective axon radius.
%   If specified, the signal scaling factor (beta) can be saved as an 
%   additional NIfTI file.
%
% EXAMPLE USAGE:
%   compute_effective_radius_multi_shell_dwi('dwi.nii.gz', 'bvals.txt', 
%       'r_eff.nii.gz', 15, 30, 2, 'Neuman', 'mask_file', 'mask.nii.gz', ...
%       'beta_file', 'beta.nii.gz', 'min_bval', 6000, 'max_bval', 30450);

    arguments
        powder_average_file string
        bval_file string
        output_file string
        delta double
        Delta double
        D0 double   
        model string
        options.mask_file (1,1) string = "" 
        options.beta_file (1,1) string = ""
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

    delta = repmat(delta, size(bvals));
    Delta = repmat(Delta, size(bvals));
    g = sqrt(bvals./(delta.^2.*(267.513*10^(-6)).^2.*(Delta-delta/3)));

    % Load mask if provided
    if options.mask_file ~= ""
        mask = niftiread(options.mask_file) > 0; % Ensure binary mask
    else
        mask = true(size(powder_average_image, 1, 2, 3)); % Use all voxels if no mask
    end
    
    % Get indices of valid mask voxels
    mask_idx = find(mask);
    M = numel(mask_idx); % Number of valid voxels
    N = size(powder_average_image, 4); % Number of values per voxel
    
    % Reshape to extract only valid voxels
    powder_average_values = reshape(powder_average_image, [], N); % Flatten first three dims
    powder_average_values = powder_average_values(mask_idx, :); % Keep only valid voxels

    % Compute r_eff_values for valid voxels
    [r_eff_values, beta_values] = estimate_r_eff_multi_shell( ...
        powder_average_values, ...
        bvals, ...
        delta, ...
        Delta, ...
        D0);

    % Create image from r_eff and beta values
    r_eff_map = nan(size(mask)); 
    beta_map = nan(size(mask)); 
    r_eff_map(mask_idx) = r_eff_values;
    beta_map(mask_idx) = beta_values;

    % write r_eff map to nifti file
    use_compression = false;
    if endsWith(output_file, ".nii.gz")
        use_compression = true;
        output_file = erase(output_file, ".nii.gz");
    elseif endsWith(output_file, ".nii")
        output_file = erase(output_file, ".nii");
    end
    nifti_info = niftiinfo(powder_average_file);
    nifti_info.ImageSize = nifti_info.ImageSize(1:ndims(r_eff_map));
    nifti_info.PixelDimensions = nifti_info.PixelDimensions(1:ndims(r_eff_map));
    nifti_info.Datatype = class(r_eff_map);
    niftiwrite(r_eff_map, ...
        output_file, ...
        nifti_info, ...
        'Compressed', ...
        use_compression);
    
    % write beta map to nifti file if desired
    if options.beta_file ~= ""
        if endsWith(options.beta_file, ".nii.gz")
            use_compression = true;
            options.beta_file = erase(options.beta_file, ".nii.gz");
        elseif endsWith(image_file, ".nii")
            options.beta_file = erase(options.beta_file, ".nii");
        end
        niftiwrite(beta_map, ...
            options.beta_file, ...
            nifti_info, ...
            'Compressed', ...
            use_compression);
    end
end

