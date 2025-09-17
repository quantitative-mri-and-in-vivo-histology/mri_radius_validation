% ========================================================================
% Illustrate histology ROI coordinates in OMM space (on T1w template).
% ========================================================================

%% compute r_eff for simulations
roi_info_table = readtable(fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/roiinfo.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");

% convert OMM coordinates to matlab's one-based indexing
roi_info_table.omm_voxel_x = roi_info_table.omm_voxel_x + 1;
roi_info_table.omm_voxel_y = roi_info_table.omm_voxel_y + 1;
roi_info_table.omm_voxel_z = roi_info_table.omm_voxel_z + 1;

% load T1w template
t1w_template_file = fullfile(getenv("MRV_DATA_PATH"), sprintf( ...
    "mri_in_vivo/templates/anat/space-omm_T1w.nii.gz"));
t1w_template = niftiread(t1w_template_file);


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
x_plot_idx = 90;
y_plot_idx = 62:178;
z_plot_idx = 58:116;

% show T1w image
nexttile(layout);
hold on;
imagesc(squeeze(t1w_template(x_plot_idx, :, :))', [0 3800]);
colormap gray;
axis image;
axis off;

% annotate ROIs on T1w image
subject_ids = unique(roi_info_table.subject_id);
for subject_index = 1:length(subject_ids)
    subject_table = roi_info_table(...
        roi_info_table.subject_id == subject_ids(subject_index),:);
    for row_index = 1:height(subject_table)
        scatter(subject_table.omm_voxel_y(row_index), ...
            subject_table.omm_voxel_z(row_index), ...
            5, ... 
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
print(gcf, '-dsvg', "rois_omm.svg");