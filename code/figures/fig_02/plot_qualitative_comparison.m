% ========================================================================
% Plot qualitative comparison of r_eff between modalities.
%
% Compute spatial patterns of r_eff across the mid-sagittal slice of the 
% corpus callosum in MNI space for histology, in-vivo and ex-vivo dMRI, 
% as well as for dMRI simulations.
% ========================================================================

%% compute r_eff for simulations
roi_info_table = readtable(fullfile(getenv("MRV_DATA_PATH"), ...
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


%% ex-vivo dMRI simulations
% load parameters and signal LUT for ex-vivo simulations
mr_protocol_ex_vivo = MrProtocol.fromJsonFile(...
    "../../../parameters/mr_protocol_ex_vivo_experimental.json");
tissue_params_ex_vivo = jsondecode(fileread( ...
    "../../../parameters/tissue_params_ex_vivo.json"));
signal_lut_ex_vivo = load(fullfile(getenv("MRV_DATA_PATH"), ...
    "/simulations/signal_per_radius_ex_vivo_experimental.mat"));

% get ex-vivo-like axon radius distributions and r_eff
[r_roi_counts_ex_vivo, r_eff_reference_ex_vivo] = get_histology_data(...
    fullfile(getenv("MRV_DATA_PATH")),...
    sim_params.r_bin_edges,...
    sim_params.radius_approximation,...
    "scaling_factor", tissue_params_ex_vivo.histology_scaling_factor);

% compute ROI signals
signals_per_roi_ex_vivo = compute_ensemble_weighted_signal( ...
    sim_params.r_bin_centers, ...
    r_roi_counts_ex_vivo, ...
    squeeze(signal_lut_ex_vivo.signal_per_r));

% simulate effective radius estimation from signals
% for experimental SNR and idealized scenario (SNR = inf)
rand_stream = RandStream(sim_params.random_generator, ...
    'Seed', sim_params.random_seed);
r_eff_simulated_idealized_ex_vivo = ...
    simulate_experimental_r_eff_estimation( ...
    signals_per_roi_ex_vivo, ...
    mr_protocol_ex_vivo, ...
    tissue_params_ex_vivo.d_0, ...
    "noise_type", "", ...
    "powder_average", "gaussian_ml", ...
    "f_im_estimate", tissue_params_ex_vivo.f_im, ...
    "rand_stream", rand_stream);
rand_stream = RandStream(sim_params.random_generator, ...
    'Seed', sim_params.random_seed);
r_eff_simulated_experiment_like_ex_vivo = ...
    simulate_experimental_r_eff_estimation( ...
    signals_per_roi_ex_vivo, ...
    mr_protocol_ex_vivo, ...
    tissue_params_ex_vivo.d_0, ...
    "noise_type", "rician", ...
    "powder_average", "rician_ml", ...
    "noise_level", 1./mr_protocol_ex_vivo.snr, ...
    "noise_level_estimate", 1./mr_protocol_ex_vivo.snr, ...
    "f_im_estimate", tissue_params_ex_vivo.f_im, ...
    "rand_stream", rand_stream);


%% in-vivo dMRI simulations
% load parameters and signal LUT for in-vivo simulations
mr_protocol_in_vivo = MrProtocol.fromJsonFile(...
    "../../../parameters/mr_protocol_in_vivo_experimental.json");
tissue_params_in_vivo = jsondecode(fileread( ...
    "../../../parameters/tissue_params_in_vivo.json"));
signal_lut_in_vivo = load(fullfile(getenv("MRV_DATA_PATH"), ...
    "/simulations/signal_per_radius_in_vivo_experimental.mat"));

% get ex-vivo-like axon radius distributions and r_eff
[r_roi_counts_in_vivo, r_eff_reference_in_vivo] = get_histology_data(...
    fullfile(getenv("MRV_DATA_PATH")),...
    sim_params.r_bin_edges,...
    sim_params.radius_approximation,...
    "scaling_factor", tissue_params_in_vivo.histology_scaling_factor);

% compute ROI signals
signals_per_roi_in_vivo = compute_ensemble_weighted_signal( ...
    sim_params.r_bin_centers, ...
    r_roi_counts_in_vivo, ...
    squeeze(signal_lut_in_vivo.signal_per_r));

% simulate effective radius estimation from signals
% for experimental SNR and idealized scenario (SNR = inf)
rand_stream = RandStream(sim_params.random_generator, ...
    'Seed', sim_params.random_seed);
r_eff_simulated_experiment_like_in_vivo = ...
    simulate_experimental_r_eff_estimation( ...
    signals_per_roi_in_vivo, ...
    mr_protocol_in_vivo, ...
    tissue_params_in_vivo.d_0, ...
    "noise_type", "rician", ...
    "powder_average", "rician_ml", ...
    "noise_level", 1./mr_protocol_in_vivo.snr, ...
    "noise_level_estimate", 1./mr_protocol_in_vivo.snr, ...
    "f_im_estimate", tissue_params_in_vivo.f_im, ...
    "rand_stream", rand_stream);
rand_stream = RandStream(sim_params.random_generator, ...
    'Seed', sim_params.random_seed);
r_eff_simulated_idealized_in_vivo = ...
    simulate_experimental_r_eff_estimation( ...
    signals_per_roi_in_vivo, ...
    mr_protocol_in_vivo, ...
    tissue_params_in_vivo.d_0, ...
    "noise_type", "", ...
    "powder_average", "gaussian_ml", ...
    "f_im_estimate", tissue_params_in_vivo.f_im, ...
    "rand_stream", rand_stream);


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
          
        % compute median r_eff within ROI mask
        mri_effective_radii = effective_radius_image(roi_mask > 0);
        ex_vivo_table.r_eff_mri_median(row_index) = ...
            median(mri_effective_radii, 'omitnan');
    end
end


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


%% generate spatial patterns across the corpus callosum for all modalities
mni_x_midslice_index = roi_info_table.mni_voxel_x(1);
in_vivo_mask = squeeze(cc_atlas_mask(mni_x_midslice_index,:,:) > 0);
plot_mni_z_idx = 83:162;
plot_mni_y_idx = 67:100;
clims_in_vivo = [1 3.5];
clims_ex_vivo = ...
    clims_in_vivo./tissue_params_in_vivo.histology_scaling_factor;
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
r_eff_map_histo_per_donor = ...
    nan([length(subject_ids), pattern_size]);
for sample_index = 1:length(subject_ids)    
    subject_table = roi_info_table(...
        roi_info_table.subject_id == subject_ids(sample_index), :);
    r_eff_map_histo_per_donor(sample_index,:,:) = ... 
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_reference_ex_vivo(subject_table.unique_id) ...
        );
end

% compute patterns per donor for in-vivo-like histology and dMRI
% simulations
r_eff_map_ex_vivo_histo_per_donor = ...
    nan([length(subject_ids), pattern_size]);
r_eff_map_ex_vivo_sim_idealized_per_donor = ...
    nan([length(subject_ids), pattern_size]);
r_eff_map_ex_vivo_sim_experiment_like_per_donor = ...
    nan([length(subject_ids), pattern_size]);

r_eff_map_in_vivo_histo_per_donor = ...
    nan([length(subject_ids), pattern_size]);
r_eff_map_in_vivo_sim_idealized_per_donor = ...
    nan([length(subject_ids), pattern_size]);
r_eff_map_in_vivo_sim_experiment_like_per_donor = ...
    nan([length(subject_ids), pattern_size]);

for sample_index = 1:length(subject_ids)    
    subject_table = roi_info_table(...
        roi_info_table.subject_id == subject_ids(sample_index), :);

    r_eff_map_ex_vivo_histo_per_donor(sample_index,:,:) = ... 
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_reference_ex_vivo(subject_table.unique_id) ...
        );
    r_eff_map_ex_vivo_sim_idealized_per_donor(sample_index,:,:) = ...
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_simulated_idealized_ex_vivo(subject_table.unique_id) ...
        );
    r_eff_map_ex_vivo_sim_experiment_like_per_donor(sample_index,:,:) = ...
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_simulated_idealized_ex_vivo(subject_table.unique_id) ...
        );

    r_eff_map_in_vivo_histo_per_donor(sample_index,:,:) = ... 
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_reference_in_vivo(subject_table.unique_id) ...
        );
    r_eff_map_in_vivo_sim_idealized_per_donor(sample_index,:,:) = ...
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_simulated_idealized_in_vivo(subject_table.unique_id) ...
        );
    r_eff_map_in_vivo_sim_experiment_like_per_donor(sample_index,:,:) = ...
        interpolate_cc( ...
            subject_table.mni_voxel_y, ...
            subject_table.mni_voxel_z, ...
            r_eff_simulated_idealized_in_vivo(subject_table.unique_id) ...
        );
end

% compute ex-vivo-like across-donor mean histology pattern
r_eff_map_ex_vivo_histo_avg = ...
    squeeze(mean(r_eff_map_ex_vivo_histo_per_donor, 1, 'omitnan'));

% extract ex-vivo patterns for ev01
r_eff_map_mri_experimental_ex_vivo = interpolate_cc( ...
    ex_vivo_table.mni_voxel_y, ...
    ex_vivo_table.mni_voxel_z,  ...
    ex_vivo_table.r_eff_mri_median(ex_vivo_table.unique_id));

% compute ex-vivo across-donor mean pattern for dMRI simulations
r_eff_map_ex_vivo_sim_idealized_avg = squeeze( ...
    mean(r_eff_map_ex_vivo_sim_idealized_per_donor, 1, 'omitnan'));
r_eff_map_ex_vivo_sim_experiment_like_avg = squeeze( ...
    mean(r_eff_map_ex_vivo_sim_experiment_like_per_donor, 1, 'omitnan'));

% compute in-vivo-like across-donor mean histology pattern
r_eff_map_in_vivo_histo_avg = ...
    squeeze(mean(r_eff_map_in_vivo_histo_per_donor, 1, 'omitnan'));

% compute in-vivo across-subject mean pattern for dMRI experiments
r_eff_map_in_vivo_experiments_avg = ...
    squeeze(r_eff_image_mean_mni_in_vivo( ...
        mni_x_midslice_index,plot_mni_z_idx,plot_mni_y_idx));

% compute in-vivo across-donor mean pattern for dMRI simulations
r_eff_map_in_vivo_sim_idealized_avg = squeeze( ...
    mean(r_eff_map_in_vivo_sim_idealized_per_donor, 1, 'omitnan'));
r_eff_map_in_vivo_sim_experiment_like_avg = squeeze( ...
    mean(r_eff_map_in_vivo_sim_experiment_like_per_donor, 1, 'omitnan'));


%% figure: histological patterns
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf, ...
    'units', 'centimeters', ...
    'position',[0, 0, 18.4, 5.0]);
t_histo = tiledlayout(1, 2, ...
    "TileSpacing", "tight", ...
    "Padding", "tight");

% plot histology patterns per subject
for subject_index = 1:length(subject_ids)
    current_ax = nexttile(t_histo);
    mean_r_eff_map_plot = ...
        squeeze(r_eff_map_histo_per_donor(subject_index,:,:))';
    plot_pattern(current_ax, ...
        mean_r_eff_map_plot, ...
        cc_mask_plot, ...
        clims_ex_vivo);
end

% add colorbar
create_cbar(current_ax, clims_ex_vivo);

% save figure
print(gcf, '-dsvg', "qualitative_comparison_histology.svg");


%% figure: ex-vivo(-like) patterns
% set up figure
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf, ...
    'units', 'centimeters', ...
    'position',[0, 0, 18.4, 3.0]);
t_ex_vivo = tiledlayout(1, 4, ...
    "TileSpacing", "none", ...
    "Padding", "tight");

% ex-vivo-like histology
current_ax = nexttile(t_ex_vivo);
plot_pattern(current_ax, ...
    r_eff_map_ex_vivo_histo_avg', ...
    cc_mask_plot, ...
    clims_ex_vivo);

% ex-vivo dMRI experiments histology
current_ax = nexttile(t_ex_vivo);
plot_pattern(current_ax, ...
    r_eff_map_mri_experimental_ex_vivo', ...
    cc_mask_plot, ...
    clims_ex_vivo);

% ex-vivo dMRI simulation (idealized)
current_ax = nexttile(t_ex_vivo);
plot_pattern(current_ax, ...
    r_eff_map_ex_vivo_sim_idealized_avg', ...
    cc_mask_plot, ...
    clims_ex_vivo);

% ex-vivo dMRI simulation (experiment-like)
current_ax = nexttile(t_ex_vivo);
plot_pattern(current_ax, ...
    r_eff_map_ex_vivo_sim_experiment_like_avg', ...
    cc_mask_plot, ...
    clims_ex_vivo);

% add colorbar
cb = create_cbar(current_ax, clims_ex_vivo);
cb.Ticks = 1:0.5:2.5;

% save figure
print(gcf, '-dsvg', "qualitative_comparison_ex_vivo.svg");


%% figure: in-vivo(-like) patterns
% set up figure
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf, ...
    'units', 'centimeters', ...
    'position',[0, 0, 18.4, 3.0]);
t_in_vivo = tiledlayout(1, 4, ...
    "TileSpacing", "none", ...
    "Padding", "tight");

% in-vivo-like histology
current_ax = nexttile(t_in_vivo);
plot_pattern(current_ax, ...
    r_eff_map_in_vivo_histo_avg', ...
    cc_mask_plot, ...
    clims_in_vivo);

% in-vivo dMRI experiments histology
current_ax = nexttile(t_in_vivo);
plot_pattern(current_ax, ...
    r_eff_map_in_vivo_experiments_avg', ...
    cc_mask_plot, ...
    clims_in_vivo);

% in-vivo dMRI simulation (idealized)
current_ax = nexttile(t_in_vivo);
plot_pattern(current_ax, ...
    r_eff_map_in_vivo_sim_idealized_avg', ...
    cc_mask_plot, ...
    clims_in_vivo);

% in-vivo dMRI simulation (experiment-like)
current_ax = nexttile(t_in_vivo);
plot_pattern(current_ax, ...
    r_eff_map_in_vivo_sim_experiment_like_avg', ...
    cc_mask_plot, ...
    clims_in_vivo);

% add colorbar
cb = create_cbar(current_ax, clims_in_vivo);
cb.Ticks = 1:0.5:3.5;

% save figure
print(gcf, '-dsvg', "qualitative_comparison_in_vivo.svg");



%% helper functions for plotting patterns
function plot_pattern(ax_handle, image, mask, clims)
    axes(ax_handle);
    hold on;
    axis image;
    axis off;
    image(~mask) = nan;
    valid_data_mask = mask & ~isnan(image);
    h = imagesc(image);
    set(h, 'AlphaData', valid_data_mask>0);
    contour(mask, [0.5 0.5], 'k', 'LineWidth', 0.5);
    clim(clims);
end

function cb = create_cbar(ax_handle, clims)
    axes(ax_handle);
    hold on;
    cb = colorbar();  % Set the colorbar position
    clim(clims);  % Set color limits based on 'clims'
    cb.Label.Interpreter = "latex";
    cb.Label.String = "$r_{\mathrm{eff}}$ [$\mu$m]";
    cb.Label.FontSize = 10;
end
