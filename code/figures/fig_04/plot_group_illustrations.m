% ========================================================================
% Illustrate axon radius distribution and group mean effective radii
% for ASD and control groups
% ========================================================================

%% set up
% load parameters and precomputed signals per axon radius table
roi_info_table = readtable( ...
    fullfile(getenv("MRV_DATA_PATH"), ...
    "histology/rawdata/roiinfo.tsv"), ...
    "FileType", "text", ...
    "Delimiter", "\t", ...
    "TextType", "string");
% extract splenium voxels
splenium_voxel_ids = find(startsWith(roi_info_table.roi_id, "S")); 
sim_params = get_default_simulation_parameters();
mr_protocol_in_vivo = MrProtocol.fromJsonFile(...
    "../../../parameters/mr_protocol_in_vivo_experimental.json");
tissue_params_in_vivo = jsondecode(fileread( ...
    "../../../parameters/tissue_params_in_vivo.json"));
signal_lut_in_vivo = load(fullfile(getenv("MRV_DATA_PATH"), ...
    "simulations/signal_per_radius_in_vivo_experimental.mat"));

% for reproducibility
rand_stream = RandStream( ...
    sim_params.random_generator, 'Seed', sim_params.random_seed);

% set up groups
group_names = ["Control", "ASD"];
group_size = 10;
num_voxels = 11;

% Default scaling factor for controls and 28.6% reduced scaling factor for 
% ASD group as reported by Wegiel et al., 2018 
% (https://doi.org/10.1186/s40478-018-0645-7) 
scaling_factors = ...
    tissue_params_in_vivo.histology_scaling_factor.*[1 (1-0.286)]; 


r_eff_simulated_mean = zeros(length(group_names), group_size);
r_section_counts_pooled = cell(1, length(group_names));
for group_index = 1:length(group_names)
    
    % simulate in-vivo like axon radius distributions
    [r_section_counts_per_scaling_factor_raw, ...
        r_effs_reference_per_simulation_raw] = ...
        get_histology_data( ...
            fullfile(getenv("MRV_DATA_PATH")), ...
            sim_params.r_bin_edges, ...
            sim_params.radius_approximation, ...
            "scaling_factor", scaling_factors(group_index));
    
    % extract roi counts for splenium regions
    r_section_counts = squeeze( ...
        r_section_counts_per_scaling_factor_raw(splenium_voxel_ids,:));

    % pool distribution for histogram illustration
    r_section_counts_pooled{group_index} = ...
        squeeze(sum(r_section_counts,1));
    
    % simulate signals per ROI
    signals_per_roi = compute_ensemble_weighted_signal( ...
        sim_params.r_bin_centers, ...
        r_section_counts, ...
        squeeze(signal_lut_in_vivo.signal_per_r));
    
    % simulate effective radii
    random_sample_ids = datasample( ...
        rand_stream, ...
        1:length(splenium_voxel_ids), ...
        num_voxels*group_size, ...
        'Replace', true);
    r_eff_simulated = simulate_experimental_r_eff_estimation( ...
        signals_per_roi(:,random_sample_ids), ...
        mr_protocol_in_vivo, ...
        tissue_params_in_vivo.d_0, ...
        "noise_type", "gaussian", ...
        "noise_level", 1./mr_protocol_in_vivo.snr, ...
        "powder_average", "gaussian_ml", ...
        "f_im_estimate", tissue_params_in_vivo.f_im, ...
        "n_noise_realizations", 1, ...
        "rand_stream", rand_stream);
    r_eff_simulated = reshape(r_eff_simulated, group_size, num_voxels); 

    % compute subject-mean across effective radii
    r_eff_simulated_mean(group_index,:) = ...
        mean(r_eff_simulated, 2, 'omitnan');
end


%% plot
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf,'units','centimeters','position',[0, 0, 11.0, 6]);
color_order = get_default_color_order();
color_order = color_order(4:5,:);
layout = tiledlayout(1,2);
layout.TileSpacing = "compact";
layout.Padding = "tight";

% illustration of group differences in axon radius distributions
nexttile(layout);
hold on;
for group_index = 1:length(group_names)
    % plot pooled splenium axon radius distribution
    histogram('BinEdges', sim_params.r_bin_edges', ...
        'BinCounts', r_section_counts_pooled{group_index}, ...
        'Normalization', 'probability', ...
        'FaceAlpha', 0.1, ...
        'EdgeColor', color_order(group_index,:), ...
        'FaceColor', color_order(group_index,:), ...
        'LineWidth', 1, ...
        'DisplayName', group_names(group_index));
    
    % annotate effective radius as vertical line
    r_eff = compute_r_eff_from_distribution( ...
        'BinCenters', sim_params.r_bin_centers, ...
        'BinCounts', r_section_counts_pooled{group_index}');
    xline(r_eff, '-', ...
        'Color', color_order(group_index,:), ...
        'LineWidth', 2, ...
        'DisplayName', sprintf("$r_{\\mathrm{eff}}$ (%s)", ...
            group_names(group_index)));
end
xlim([0 3.4]);
ylabel("probability", 'Interpreter', 'latex');
xlabel("$r [\mu$m$]$", 'Interpreter', 'latex');
pbaspect([1 1 1]);
box on;
lgd = legend('Interpreter', 'latex', 'NumColumns', 2, 'FontSize', 7);
lgd.Location = "northoutside";
lgd.ItemTokenSize(1) = 15;


% illustration of group differences in subject-level mean effective radii
% of the splenium
nexttile(layout);
hold on;
for group_index = 1:length(group_names)
    boxchart(repmat(categorical(group_names(group_index)), ...
        length(r_eff_simulated_mean(group_index,:)), 1), ...
        r_eff_simulated_mean(group_index,:), ['' ...
        'BoxFaceColor'], color_order(group_index,:));
end
box on;
xlabel('group', 'Interpreter', 'latex');
ylabel('$\bar{r}_{\mathrm{eff}}$ [$\mu$m]', 'Interpreter', 'latex');
pbaspect([1 1 1]);

% add signifcance marker
x1 = 1; % x-coordinate of the first box
x2 = 2; % x-coordinate of the second box
y1 = max(r_eff_simulated_mean(1,:)); % Top position of the first box
y2 = max(r_eff_simulated_mean(2,:)); % Top position of the second box
y_max = max(y1, y2) + 0.1; % Position for the line and asterisk 
line([x1, x1, x2, x2], ...
    [y1 + 0.05, y_max, y_max, y2 + 0.05], ...
    'Color', 'k', ...
    'HandleVisibility', 'off');
text((x1 + x2) / 2, y_max + 0.05, ...
    '$p < 0.05$?', ...
    'HorizontalAlignment', 'center', ...
    'FontSize', 8, ...
    'Interpreter', 'latex');
ylims = ylim();
ylims(1) = ylims(1) - 0.05; 
ylims(2) =  y_max + 0.15;
ylim(ylims);

% save figure
print(gcf, '-dsvg', "group_illustrations.svg");
