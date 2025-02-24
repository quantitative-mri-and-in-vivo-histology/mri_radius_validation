% ========================================================================
% Plot sampling statistics of our dataset versus existing corpus callosum
% datasets.
% ========================================================================

% load data
r_bin_centers = readtable(fullfile(getenv("MRV_DATA_PATH"), ... 
    "histology/rawdata/desc-binCenters_radii.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t");
r_bin_centers = r_bin_centers.Variables;
r_bin_edges = readtable(fullfile(getenv("MRV_DATA_PATH"), ... 
    "histology/rawdata/desc-binEdges_radii.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t");
r_bin_edges = r_bin_edges.Variables;
r_bin_counts = readtable(fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/desc-countsCircularEq_radii.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t");
r_bin_counts = r_bin_counts.Variables;
roi_info_table = readtable(fullfile(getenv("MRV_DATA_PATH"), ... 
    "histology/rawdata/roiinfo.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");

% prepare statistics
subject_ids = unique(roi_info_table.subject_id);
samples_per_donor_ours = zeros(length(subject_ids),1);
for donor_index = 1:length(subject_ids)
    samples_per_donor_ours(donor_index) = ...
        sum(roi_info_table.subject_id == subject_ids(donor_index));
end
study_names = [
    "Aboiitz et al.";
    "Caminiti et al."; 
    "Liewald et al.";
    "Ours"];
donors_per_dataset = [1, 3, 2, 2];
rois_per_donors = [5, 4, 3, mean(samples_per_donor_ours)];
areas_per_roi = [45*45*10^-6, 112*87*10^(-6), 13*13*8*10^(-6), 8];
axons_per_roi = [300, 1750, 300, 10^6];

% set up subplots
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf,'units','centimeters','position',[0, 0, 12.1, 5.3]);
layout = tiledlayout(1, 3, ...
    "TileSpacing", "loose", ...
    "Padding", "tight");
color_order = get_default_color_order();
color_order = [color_order; 0,0.47,0.84];
plot_axes = gobjects(1,3);
plot_axes(1) = nexttile(layout);
plot_axes(2) = nexttile(layout);
plot_axes(3) = nexttile(layout);

for study_index = 1:length(study_names)
   
    % plot donors per dataset
    axes(plot_axes(1));
    hold on;
    bar(categorical(study_names(study_index)), ...
        donors_per_dataset(study_index), ...
        'FaceColor', color_order(study_index,:), ...
        'LineWidth', 1,...
        'DisplayName', study_names(study_index));
    ylabel("donors", 'interpreter', 'latex');
    yticks([0:1:max(donors_per_dataset)]);
    xlabel("dataset", 'Interpreter', 'latex');
    ylim([0, max(donors_per_dataset)+0.1]);
    xticklabels([]);
    pbaspect([1 1 1]);
    box on;
    
    % plot rois per donor
    axes(plot_axes(2));
    hold on;
    bar(categorical(study_names(study_index)), ...
        rois_per_donors(study_index), ...
        'FaceColor', color_order(study_index,:), ...
        'LineWidth', 1);
    ylabel("ROIs per donor", 'Interpreter', 'latex');
    xlabel("dataset", 'Interpreter', 'latex');
    ylim([0, max(rois_per_donors)+0.5]);
    xticklabels([]);
    pbaspect([1 1 1]);
    box on;

    % plot area/axon count per roi
    axes(plot_axes(3));
    hold on;
    xlim([0 1.5*10^1]);
    ylim([0 5*10^6]);
    xline(5.25, 'k--', 'LineWidth', 2);
    plot(areas_per_roi(study_index), axons_per_roi(study_index), ...
        'o', ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', color_order(study_index,:), ...
        'MarkerSize', 7);
    xlabel("area per ROI [$\mu$m$^2$]", 'Interpreter', 'latex');
    ylabel("axons per ROI", 'Interpreter', 'latex');
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    yticks(logspace(1, 8, 8));
    xticks(logspace(-4, 1, 6));
    pbaspect([1 1 1]);
    box on;
end

% add global dataset legend
axes(plot_axes(1));
lgd = legend('NumColumns', 4, ...
    'Interpreter', 'latex');
lgd.Layout.Tile = "North";

% save figure
print(gcf, '-dsvg', "dataset_comparison.svg");