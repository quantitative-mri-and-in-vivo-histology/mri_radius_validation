% ========================================================================
% Plot axon radius distribution comparison between our data
% and a subsample mimicking the sample size of Aboitiz et al.
% (see https://doi.org/10.1016/0006-8993(92)90178-C)
% ========================================================================

% load data
r_bin_centers = readtable(fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/desc-binCenters_radii.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t");
r_bin_centers = r_bin_centers.Variables;
r_bin_counts = readtable(fullfile(getenv("MRV_DATA_PATH"), ... 
    "histology/rawdata/desc-countsCircularEq_radii.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t");
r_bin_counts = r_bin_counts.Variables;
sim_params = get_default_simulation_parameters();

% extract radii of exemplary ROI
example_section_index = 33;
radii_full_sample = convert_bin_counts_to_set( ...
    r_bin_counts(example_section_index, :), r_bin_centers);

% subsample radii to generate Aboitiz-like subsample of 1000 axons
% (as in https://doi.org/10.1016/0006-8993(92)90178-C)
aboitiz_sample_size = 10^3;
rng(481, sim_params.random_generator);
radii_aboitiz_like = datasample(radii_full_sample, aboitiz_sample_size, ...
    'Replace', false);
radii = {radii_full_sample, radii_aboitiz_like};

% set up plots
fig_handle = figure;
layout = tiledlayout(1, 2, ...
    'Padding', 'tight', ...
    'TileSpacing', 'tight');
set(gcf, ...
    'units', 'centimeters', ...
    'position', [0, 0, 8.4, 5]);
plot_axes = gobjects(1, 2);
lm_color = [0, 0.47, 0.84];
color_order = get_default_color_order();
em_color = color_order(1,:);
colors = [lm_color; em_color];

% plot axon radius distributions
for plot_index = 1:length(radii)
    plot_axes(plot_index) = nexttile(layout);
    hold on;
    histogram(radii{plot_index}, ...
        0:0.1:5, ...
        'Normalization', 'pdf', ...
        'LineWidth', 0.7, ...
        'FaceColor', colors(plot_index,:), ...
        'FaceAlpha', 1, ...
        'EdgeColor', [0 0 0]);
    xline(compute_r_eff_from_distribution(radii{plot_index}), ':', ...
        'Color', colors(plot_index,:), ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off');
    xlabel("$r~[\mu$m$]$", 'Interpreter', 'latex');
    ylabel("probability", 'Interpreter', 'latex');
    pbaspect([1 1 1]);
    xlims = [0 4.1];
    xlim(xlims);
    xticks(xlims(1):1:xlims(2));
    ylim([0 2.2]);
    box on;
end

% plot insets illustrating the tail of the distribution
for plot_index = 1:length(radii)
    inset_xlims = [1.5 4.1];
    inset_ylims = [0.0 0.085];
    inset_ax = plot_inset(plot_axes(plot_index), ...
        [0.40 0.43 0.58 0.58], ...
        inset_xlims, ...
        inset_ylims);
    yticks(inset_ax, 0:0.04:0.08);
    axes(plot_axes(plot_index));
end

% save figure
print(gcf, '-dsvg', "histogram_examples.svg");