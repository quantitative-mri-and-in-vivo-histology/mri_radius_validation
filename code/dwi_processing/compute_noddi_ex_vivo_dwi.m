function compute_noddi_ex_vivo_dwi( ...
    image_file, ...
    bvec_file, ...
    bval_file, ...
    outfile_prefix, ...
    noisemap_file, ...
    diffusivity, ...
    options)
% Estimates NODDI (Neurite Orientation Dispersion and Density Imaging) 
% parameters from diffusion-weighted MRI data. 
%
% USAGE:
%   compute_noddi_dwi(image_file, bvec_file, bval_file, outfile_prefix, ...
%       noisemap_file, "exvivo", 'mask_file', "brain_mask.nii.gz", ...
%       'diffusivity', 3.5e-10)
%
% INPUTS:
%   image_file     - (string) Path to the input diffusion-weighted image 
%                    (4D NIfTI file).
%   bvec_file      - (string) Path to the corresponding b-vector file.
%   bval_file      - (string) Path to the corresponding b-value file.
%   outfile_prefix - (string) Prefix for output files (including directory).
%   noisemap_file  - (string) Path to a noise map NIfTI file.
%   model          - (string) NODDI model type: 
%                    - "invivo"  -> WatsonSHStickTortIsoV_B0
%                    - "exvivo"  -> WatsonSHStickTortIsoVIsoDot_B0
%
% OPTIONS (Name-Value Pairs):
%   options.diffusivity - (numeric, default = []) If specified, 
%                         constrains the diffusivity parameter 
%                         in the NODDI model.
%   options.mask_file   - (string, default = "") Path to a binary mask 
%                         NIfTI file. If not provided, a dummy mask 
%                         is created and saved.
%
% OUTPUTS:
%   - The function generates the following output files:
%     1. '<outfile_prefix>_roi.mat' - ROI data for fitting.
%     2. '<outfile_prefix>_fitted_params.mat' - Fitted NODDI parameters.
%     3. NIfTI files containing the final computed NODDI parameter maps.

    arguments
        image_file (1,1) string
        bvec_file (1,1) string
        bval_file (1,1) string
        outfile_prefix (1,1) string
        noisemap_file (1,1) string
        diffusivity (1,1) double
        options.mask_file (1,1) string = "" 
    end
    
    noisemap = niftiread(noisemap_file);

    if options.mask_file == ""
        % create and save dummy mask
        mask_image = ones(size(noisemap));
        options.mask_file_stripped = sprintf("%s_mask", outfile_prefix);
        options.mask_file = sprintf("%s.nii", options.mask_file_stripped);
        nifti_info = niftiinfo(noisemap_file);
        nifti_info.ImageSize = nifti_info.ImageSize(1:ndims(mask_image));
        nifti_info.PixelDimensions = ...
            nifti_info.PixelDimensions(1:ndims(mask_image));
        nifti_info.Datatype = class(mask_image);
        niftiwrite(mask_image, options.mask_file_stripped, ...
            nifti_info, 'Compressed', false);     
    end

    image_file = decompress_if_needed(image_file, outfile_prefix);
    if options.mask_file ~= ""
        options.mask_file = decompress_if_needed(options.mask_file, ...
            outfile_prefix);
    end

    noddi_roi_file = sprintf("%s_roi.mat", outfile_prefix);
    fitted_params_file = sprintf("%s_fitted_params.mat", outfile_prefix);

    CreateROI(char(image_file), ...
        char(options.mask_file ), ...
        char(noddi_roi_file));
    protocol = FSL2Protocol(char(bval_file), char(bvec_file));

    noddi_model = MakeModel('WatsonSHStickTortIsoVIsoDot_B0');
    mask_image = niftiread(options.mask_file);
    global_sigma = double(mean(noisemap(mask_image>0), 'all'));
    noddi_model.sigma.perVoxel=false;
    noddi_model.sigma.globalSigma=global_sigma;

    noddi_model.GS.fixedvals(2) = diffusivity;
    noddi_model.GD.fixedvals(2) = diffusivity;
    
    batch_fitting( ...
        char(noddi_roi_file), ...
        protocol, ...
        noddi_model, ...
        char(fitted_params_file));
    SaveParamsAsNIfTI( ...
        char(fitted_params_file), ...
        char(noddi_roi_file), ...
        char(options.mask_file), ...
        char(outfile_prefix));
end


function uncompressed_file = decompress_if_needed(file, outfile_prefix)
    % Check if the file ends with .nii.gz and decompress it
    out_dir = fileparts(outfile_prefix);
    if endsWith(file, ".nii.gz")
        uncompressed_files = gunzip(file, out_dir);
        uncompressed_file = uncompressed_files{1};
    else
        uncompressed_file = file;
    end
end