% ========================================================================
% Illustrate histology ROI coordinates in native in-vivo dMRI space 
% (on T1w image of exemplary in-vivo subject).
% ========================================================================

%% compute r_eff for simulations
roi_info_table = readtable(fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/roiinfo.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");
sim_params = get_default_simulation_parameters();


%% in-vivo dMRI experiments
% read T1w image and transform to voxel coordinates
subject_id_index = 1;
in_vivo_processed_data_dir = fullfile(getenv("MRV_DATA_PATH"), ...
    "mri_in_vivo/processed");
processed_subj_dir_fstructs = dir( ...
    fullfile(in_vivo_processed_data_dir, "sub*"));
subject_id = sscanf( ...
    processed_subj_dir_fstructs(subject_id_index).name, 'sub-%s');
t1w_image_file = fullfile(in_vivo_processed_data_dir, ...
    sprintf("sub-%s/anat/sub-%s_T1w.nii.gz", ...
    subject_id, subject_id));
t1w_image = niftiread(t1w_image_file);
t1w_info = niftiinfo(t1w_image_file);
t1w_T = inv(t1w_info.Transform.T');

% read physical coordinates of histological ROIs in native MRI space
loc_table = readtable(fullfile(in_vivo_processed_data_dir, ...
    sprintf("sub-%s/dwi/sub-%s_roiinfo.tsv", ...
    subject_id, subject_id)), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");
coord_table = loc_table(:, ...
    ["native_physical_x", "native_physical_y", "native_physical_z"]);  
coord_physical = coord_table.Variables;

% compute T1w voxel coordinates from physical coordinates (ITK conv.)
coord_physical(:,1:2) = coord_physical(:,1:2).*(-1);
coord_voxel = t1w_T(1:3,1:3)*coord_physical' + t1w_T(1:3,4);
coord_voxel = int16(round(coord_voxel + 1));


%% figure
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf, ...
    'units', 'centimeters', ...
    'position',[0, 0, 7, 5.0]);
layout = tiledlayout(1, 1, ...
    "TileSpacing", "tight", ...
    "Padding", "tight");
color_order = get_default_color_order();
y_plot_idx = 60:210;
z_plot_idx = 130:230;
x_plot_idx = 96;

% show T1w image
nexttile(layout);
hold on;
imagesc(squeeze(t1w_image(x_plot_idx, :, :))', [0 300]);
colormap gray;
axis image;
axis off;

% annotate ROIs on T1w image
subject_ids = unique(loc_table.subject_id);
for subject_index = 1:length(subject_ids)
    coord_table = loc_table(...
        loc_table.subject_id == subject_ids(subject_index),:);
    for row_index = 1:height(coord_table)
        scatter(coord_voxel(2, coord_table.unique_id(row_index)), ...
            coord_voxel(3, coord_table.unique_id(row_index)), ...
            3, ... 
            color_order(subject_index,:), ... 
            'filled', ... 
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.7);
    end
end

% limit view to corpus callosum ROI
xlim(minmax(y_plot_idx));
ylim(minmax(z_plot_idx));

% save_figure
print(gcf, '-dsvg', "rois_native.svg");