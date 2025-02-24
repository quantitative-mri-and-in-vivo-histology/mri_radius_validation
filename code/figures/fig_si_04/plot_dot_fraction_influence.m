% ========================================================================
% Assess influence of dot fraction estimates on ex-vivo dMRI r_eff
% ========================================================================

%% compute r_eff for simulations
roi_info_table = readtable( ...
    fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/roiinfo.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");
sim_params = get_default_simulation_parameters();
cc_eval_params = jsondecode(fileread(...
    "../../../parameters/cc_evaluation_parameters.json"));

% convert MNI coordinates to matlab's one-based indexing
roi_info_table.mni_voxel_x = roi_info_table.mni_voxel_x + 1;
roi_info_table.mni_voxel_y = roi_info_table.mni_voxel_y + 1;
roi_info_table.mni_voxel_z = roi_info_table.mni_voxel_z + 1;

%% in-vivo dMRI experiments
subject_ids = unique(roi_info_table.subject_id);
cc_atlas_mask_file = fullfile(getenv("MRV_DATA_PATH"), ...
    "mri_in_vivo/templates/masks", ...
    "space-MNI152NLin6Asym_label-ccMid_mask.nii.gz");
cc_atlas_mask = niftiread(cc_atlas_mask_file);

in_vivo_processed_data_dir = fullfile(getenv("MRV_DATA_PATH"), ...
    "mri_in_vivo/processed");
processed_subj_dir_fstructs = dir( ...
    fullfile(in_vivo_processed_data_dir, "sub*"));

in_vivo_r_eff_images = ...
    nan([length(processed_subj_dir_fstructs), size(cc_atlas_mask)]);

% read per-subject effective radius map in MNI space and segment CC
for subject_id_index = 1:length(processed_subj_dir_fstructs)
    subject_id = sscanf( ...
        processed_subj_dir_fstructs(subject_id_index).name, 'sub-%s');

    % read subject r_eff and FA images
    r_eff_image = niftiread(fullfile(in_vivo_processed_data_dir, ...
        sprintf("sub-%s/dwi/" + ...
        "sub-%s_space-MNI152NLin6Asym_effectiveRadius.nii.gz", ...
        subject_id, subject_id)));
    fa_image =  niftiread(fullfile(in_vivo_processed_data_dir, ...
        sprintf("sub-%s/dwi/sub-%s_space-MNI152NLin6Asym_FA.nii.gz", ...
        subject_id, subject_id)));

    % exclude values outside corpus callosum
    cc_mask = cc_atlas_mask > 0 ...
        & fa_image >= cc_eval_params.fa_min ...
        & r_eff_image >= cc_eval_params.r_eff_min;
    r_eff_image(~cc_mask) = nan;
    in_vivo_r_eff_images(subject_id_index,:,:,:) = r_eff_image;
end

n_valid_subjects = squeeze(sum(~isnan(in_vivo_r_eff_images),1));
valid_voxels_mask = n_valid_subjects >= 3;
val_idx_rep = repmat(shiftdim(valid_voxels_mask, -1), ...
    size(in_vivo_r_eff_images,1), 1, 1, 1);
in_vivo_r_eff_images(~val_idx_rep) = nan;
r_eff_image_mean_mni_in_vivo = squeeze( ...
    mean(in_vivo_r_eff_images, 1, 'omitnan'));


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
        roi_id = sscanf(matches{3}, '%s');
        
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



%% generate spatial patterns across the corpus callosum for all modalities
cc_mask_image_file = fullfile(getenv("MRV_DATA_PATH"), ...
    "mri_in_vivo/templates/masks", ...
    "space-MNI152NLin6Asym_label-ccMid_mask.nii.gz");
cc_mask_image = niftiread(cc_mask_image_file);
mni_x_midslice_index = roi_info_table.mni_voxel_x(1);
in_vivo_mask = squeeze(cc_mask_image(mni_x_midslice_index,:,:) > 0);
plot_mni_z_idx = 83:162;
plot_mni_y_idx = 67:100;
cc_mask_plot = in_vivo_mask ...
    & squeeze(valid_voxels_mask(mni_x_midslice_index, :, :));
cc_mask_plot = cc_mask_plot(plot_mni_z_idx,plot_mni_y_idx)';

% interpolation function to generate corpus callosum patterns
[xGrid, yGrid] = meshgrid(plot_mni_z_idx, ...
    plot_mni_y_idx);
interpolate_cc = @(x, y, vals) ...
    griddata(x, y, vals, xGrid, yGrid, 'nearest')'; 

% compute histological patterns per donor
pattern_size = [length(plot_mni_z_idx), length(plot_mni_y_idx)];

% extract ex-vivo patterns for ev01
r_eff_map_mri_experimental_ex_vivo = interpolate_cc( ...
    ex_vivo_table.mni_voxel_y, ...
    ex_vivo_table.mni_voxel_z,  ...
    ex_vivo_table.r_eff_mri_median(ex_vivo_table.unique_id));

% extract ex-vivo patterns for ev01
dot_fraction_mri_experimental_ex_vivo = interpolate_cc( ...
    ex_vivo_table.mni_voxel_y, ...
    ex_vivo_table.mni_voxel_z,  ...
    ex_vivo_table.dot_fraction_mri_median(ex_vivo_table.unique_id));


%% figure: r_eff and dot fraction patterns
% set up figure
color_order = get_default_color_order();
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf, ...
    'units', 'centimeters', ...
    'position',[0, 0, 6 5.5]);
layout_x_y_plot = tiledlayout(2,1, ...
    "TileSpacing", "tight", ...
    "Padding", "none");

% plot effective radius pattern
clims_r_eff = [0.5 1.5];
nexttile(layout_x_y_plot);
hold on;
axis image;
axis off;
image = r_eff_map_mri_experimental_ex_vivo';
image(~cc_mask_plot) = nan;
valid_data_mask = cc_mask_plot & ~isnan(image);
h = imagesc(image);
set(h, 'AlphaData', valid_data_mask>0);
contour(cc_mask_plot, [0.5 0.5], 'k', 'LineWidth', 0.5);
cb = colorbar();  % Set the colorbar position
clim(clims_r_eff);  % Set color limits based on 'clims'
cb.Label.Interpreter = "latex";
cb.Label.String = "$r_{\mathrm{eff}}$ [$\mu$m]";
cb.Label.FontSize = 8;

% plot dot fraction pattern
clims_dot = [0.1 0.4];
nexttile(layout_x_y_plot);
hold on;
axis image;
axis off;
image = dot_fraction_mri_experimental_ex_vivo';
image(~cc_mask_plot) = nan;
valid_data_mask = cc_mask_plot & ~isnan(image);
h = imagesc(image);
set(h, 'AlphaData', valid_data_mask>0);
contour(cc_mask_plot, [0.5 0.5], 'k', 'LineWidth', 0.5);
cb = colorbar();  % Set the colorbar position
clim(clims_dot);  % Set color limits based on 'clims'
cb.Ticks = min(clims_dot):0.1:max(clims_dot);
cb.Label.Interpreter = "latex";
cb.Label.String = "$f_{\mathrm{im}}$";
cb.Label.FontSize = 8;

% save figure
print(gcf, '-dsvg', "dot_fraction_influence_patterns.svg");


%% figure: dot fraction versus r_eff plot
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf, ...
    'units', 'centimeters', ...
    'position',[0, 0, 5, 5.5]);
layout_x_y_plot = tiledlayout(1, 1, ...
    "TileSpacing", "tight", ...
    "Padding", "none");

region_names = ["genu", ...
    "ant. midbody", ... 
    "midbody", ...
    "post. midbody", ... 
    "splenium"];
region_prefixes = ["G", "AM", "M", "PM", "S"];

% plot r_eff versus dot fraction per CC subregion
current_ax = nexttile(layout_x_y_plot);
hold on;
for region_index = 1:length(region_names)
    
    region_idx = startsWith(ex_vivo_table.roi_id, ...
        region_prefixes(region_index));

    plot(ex_vivo_table.dot_fraction_mri_median(region_idx), ...
        ex_vivo_table.r_eff_mri_median(region_idx), ...
        'o', ...
        'MarkerSize', 5, ...
        'MarkerFaceColor', color_order(region_index,:), ...
        'MarkerEdgeColor', 'k', ...
        'DisplayName', region_names(region_index));
end
box on;
xlabel("$f_{\mathrm{im}}$", ...
    'Interpreter', 'latex');
ylabel("$r_{\mathrm{eff}}$ [$\mu$m]", ...
    'Interpreter', 'latex');
pbaspect([1 1 1]);
lgd = legend('NumColumns', 2, ...
    'Interpreter', 'latex', ...
    'FontSize', 6);
lgd.ItemTokenSize(1) = 12;
lgd.Location = "northoutside";

% save figure
print(gcf, '-dsvg', "dot_fraction_influence_quantitative.svg");

