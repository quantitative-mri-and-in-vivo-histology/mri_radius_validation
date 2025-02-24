% ========================================================================
% Compute mean dot fraction for ex-vivo dMRI data
% ========================================================================

%% compute r_eff for simulations
roi_info_table = readtable(fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/roiinfo.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");
sim_params = get_default_simulation_parameters();

% convert MNI coordinates to matlab's one-based indexing
roi_info_table.mni_voxel_x = roi_info_table.mni_voxel_x + 1;
roi_info_table.mni_voxel_y = roi_info_table.mni_voxel_y + 1;
roi_info_table.mni_voxel_z = roi_info_table.mni_voxel_z + 1;


%% ex-vivo dMRI simulations
% load parameters and signal LUT for ex-vivo simulations
mr_protocol_ex_vivo = MrProtocol.fromJsonFile( ...
    "../../parameters/mr_protocol_ex_vivo_experimental.json");
tissue_params_ex_vivo = jsondecode(fileread( ...
    "../../parameters/tissue_params_ex_vivo.json"));


%% ex-vivo dMRI experiments
% set data paths for experimental ex-vivo dMRI data
ex_vivo_data_dir = fullfile(getenv("MRV_DATA_PATH"), "mri_ex_vivo");
ex_vivo_processed_data_dir = fullfile(ex_vivo_data_dir, "processed");
ex_vivo_masks_data_dir = fullfile(ex_vivo_data_dir, "masks");

% copy ROI info table for ex-vivo data (only sub-ev01) and
% add column for r_eff from dMRI experiments
ex_vivo_table = roi_info_table(roi_info_table.subject_id == "ev01", :);
ex_vivo_table.r_eff_mri_median = ...
    nan(height(ex_vivo_table), 1);
ex_vivo_table.dot_fraction_mri_median = ...
    nan(height(ex_vivo_table), 1);

% collect r_eff from ex-vivo dMRI for each histological ROI
for row_index = 1:height(ex_vivo_table)
     
    subject_id = ex_vivo_table.subject_id(row_index);      
    roi_mask_fstructs = dir(fullfile(ex_vivo_masks_data_dir, ...
        sprintf("**/dwi/sub-%s_*roi-%s_*.nii.gz", ...
        subject_id, ...
        ex_vivo_table.roi_id(row_index))));

    if ~isempty(roi_mask_fstructs)
        % extract sub, sample and roi tags
        pattern = "sub-([^_]+)_sample-(\d+)_roi-([^_]+)";
        matches = regexp(roi_mask_fstructs(1).name, pattern, ...
            'tokens', 'once');
        sample_id = sscanf(matches{2}, '%d');
        
        % read ROI mask
        roi_mask_file = fullfile(roi_mask_fstructs(1).folder, ...
            roi_mask_fstructs(1).name);
        roi_mask = niftiread(roi_mask_file);
    
        % read effective radius image
        effective_radius_image_file = ...
            fullfile(ex_vivo_processed_data_dir, ...
            sprintf("sub-%s", subject_id), ...
            "dwi", ...
            sprintf("sub-%s_sample-%02d_effectiveRadius.nii.gz", ...
            subject_id, sample_id));
        effective_radius_image = niftiread(effective_radius_image_file);

        dot_fraction_image_file = ...
            fullfile(ex_vivo_processed_data_dir, ...
            sprintf("sub-%s", subject_id), ...
            "dwi", ...
            sprintf("sub-%s_sample-%02d_dotFraction.nii.gz", ...
            subject_id, sample_id));
        dot_fraction_image = niftiread(dot_fraction_image_file);
          
        % compute median r_eff within ROI mask
        mri_effective_radii = effective_radius_image(roi_mask > 0);
        ex_vivo_table.r_eff_mri_median(row_index) = ...
            median(mri_effective_radii, 'omitnan');

        mri_dot_fraction = dot_fraction_image(roi_mask > 0);
        ex_vivo_table.dot_fraction_mri_median(row_index) = ...
            median(mri_dot_fraction, 'omitnan');
    end
end



disp(['dot fraction = ', num2str( ...
    mean(ex_vivo_table.dot_fraction_mri_median, 'omitnan'))]);
