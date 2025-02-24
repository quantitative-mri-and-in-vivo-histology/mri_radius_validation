function compute_dot_fraction_dwi(dwi_image_file,...
    bvec_file,...
    bval_file,...
    fibre_dirs_image_file,...
    dot_fraction_image_file,...
    noisemap_file,...
    max_angle,...
    min_bval,...
    options)
% Estimates the dot fraction from diffusion-weighted MRI (dMRI) data.
%
%   This function computes the dot fraction by fitting a Rician distribution to 
%   dMRI signals highly aligned with fiber orientations. It processes only signals 
%   with a b-value ≥ `min_bval`and an angular deviation ≤ `max_angle` from 
%   the local fiber direction.
%
% INPUTS:
%   dwi_image_file         - (string) Path to the 4D NIfTI dMRI file.
%   bvec_file              - (string) Path to the gradient directions file (BVEC).
%   bval_file              - (string) Path to the gradient strengths file (BVAL).
%   fibre_dirs_image_file  - (string) Path to the NIfTI file with fiber orientations.
%   dot_fraction_image_file- (string) Path to save the computed dot fraction image.
%   noisemap_file          - (string) Path to the NIfTI noise estimate file.
%   max_angle              - (numeric) Maximum angular deviation (degrees) for inclusion.
%   min_bval               - (numeric) Minimum b-value threshold for signal selection.
%
% OPTIONAL:
%   options.mask_file      - (string) Path to a binary mask (NIfTI). If omitted, 
%                            computations are performed on the full image.
%
% OUTPUT:
%   A NIfTI file containing voxel-wise dot fraction values.
%
% EXAMPLE USAGE:
%   compute_dot_fraction_dwi('dwi.nii.gz', 'bvecs.txt', 'bvals.txt', ...
%       'fibre_dirs.nii.gz', 'dot_fraction.nii.gz', 'noise_map.nii.gz', ...
%       20, 1000, struct('mask_file', 'mask.nii.gz'));
%

    arguments
        dwi_image_file string
        bvec_file string
        bval_file string
        fibre_dirs_image_file string
        dot_fraction_image_file string
        noisemap_file string
        max_angle {mustBeNumeric}
        min_bval {mustBeNumeric}
        options.mask_file (1,1) string = ""
    end
    
    angle = @(a,b) atan2(norm(cross(a,b)), dot(a,b));
    dwi_image = niftiread(dwi_image_file);
    bvec = importdata(bvec_file)';
    bval = importdata(bval_file)';

    fiber_dirs_image = niftiread(fibre_dirs_image_file);
    noisemap_image = niftiread(noisemap_file);

    bvec = bvec(bval >= min_bval,:);
    dwi_image = dwi_image(:,:,:,bval >= min_bval);
    bval = bval(bval >= min_bval);
    bval_unique = unique(bval);

    if options.mask_file ~= ""
        mask_image = niftiread(options.mask_file);
    else
        mask_image = ones(size(dwi_image, 1:3));
    end

    mask_indices = find(mask_image > 0);
    opts = statset('mlecustom');
    opts.FunValCheck = 'off';
    dot_fraction_image = NaN(size(mask_image));
    dot_fraction_values = size(mask_indices);
    
    parfor i = 1:length(mask_indices)
        [i1,i2,i3] = ind2sub(size(mask_image), mask_indices(i));
        dot_frac_signals = [];

        % collect signals for which fiber direction aligns sufficiently
        % with main fiber direction (angle <= max_angle)
        % and b >= b_min
        for shell_index = 1:length(bval_unique)
    
            signal_values = squeeze(dwi_image(i1,i2,i3, ...
                bval == bval_unique(shell_index)));
            fibre_direction = squeeze(fiber_dirs_image(i1,i2,i3,:));
            bvec_shell = bvec(bval == bval_unique(shell_index),:);
            angles = zeros(length(bvec_shell),1);
    
            for j = 1:length(bvec_shell)
                angles(j) = rad2deg(angle(fibre_direction, ...
                    bvec_shell(j,:)'));
                if angles(j) > 90
                    angles(j) = 180-angles(j);
                end
            end
            
            if bval_unique(shell_index) >= min_bval
                dot_frac_signals = [dot_frac_signals; 
                    signal_values(angles <= max_angle)];
            end
           
        end 
        
        noise_estimate = noisemap_image(mask_indices(i));
        % fit Rician to collected signals
        rician_pdf = @(x,mu,sigma) pdf('rician', x, mu, noise_estimate);
        [phat, pci] = mle(dot_frac_signals, ...
            'pdf', rician_pdf, ...
            'start', 0.1, ...
            'Options', opts);
        dist = makedist("Rician", phat, noise_estimate);
        dot_fraction_values(i) = dist.mean;
    end
    dot_fraction_image(mask_indices) = dot_fraction_values;

    % write dot fraction map to nifti file
    use_compression = false;
    if endsWith(dot_fraction_image_file, ".nii.gz")
        use_compression = true;
        dot_fraction_image_file = ...
            erase(dot_fraction_image_file, ".nii.gz");
    elseif endsWith(image_file, ".nii")
        dot_fraction_image_file = erase(dot_fraction_image_file, ".nii");
    end
    nifti_info = niftiinfo(dwi_image_file);
    nifti_info.ImageSize = ...
        nifti_info.ImageSize(1:ndims(dot_fraction_image));
    nifti_info.PixelDimensions = ...
        nifti_info.PixelDimensions(1:ndims(dot_fraction_image));
    nifti_info.Datatype = class(dot_fraction_image);
    niftiwrite(dot_fraction_image, ...
        dot_fraction_image_file, ...
        nifti_info, 'Compressed', ...
        use_compression);
end