% ========================================================================
% Illustrate axon radius distributions of different axon populations. 
% To this end, extropolate human corpus callosum radii to the rat corpus 
% callosum and the human cortical spinal tract.
% ========================================================================

%% simulate data
% set up parameters
sim_params = get_default_simulation_parameters();
scaling_factors = [0.5 1 1.15]*1.0;
r_bin_edges = sim_params.r_bin_edges;


%% run simulations
axon_radii = cell(1, length(scaling_factors));
r_eff = zeros(1, length(scaling_factors));

for scaling_factor_index = 1:size(scaling_factors,2)
    % get scaled axon radius distributions
    r_counts = get_histology_data(...
        fullfile(getenv("MRV_DATA_PATH")),...
        sim_params.r_bin_edges,...
        sim_params.radius_approximation,...
        "scaling_factor", scaling_factors(scaling_factor_index));
    r_counts = sum(r_counts, 1);
    axon_radii{scaling_factor_index} = convert_bin_counts_to_set( ...
        r_counts, sim_params.r_bin_centers);
    r_eff(scaling_factor_index) = ...
        compute_r_eff_from_distribution(axon_radii{scaling_factor_index});
end


%% plot results
% set up subplots
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf,'units','centimeters','position',[0, 0, 5.4, 5]);
layout = tiledlayout(1, 1, ...
    "TileSpacing", 'tight', ...
    'Padding', 'tight');
color_order = get_default_color_order();
face_alphas = [0.2 0.0 0];
edge_alphas = [1 0.5 0];
nexttile(layout);
hold on;

% plot axon radius distributions
for scaling_factor_index = 1:length(scaling_factors)
    histogram(axon_radii{scaling_factor_index}, ...
        'BinEdges', r_bin_edges, ...
        'Normalization', 'probability',...
        'LineWidth', 1, ...
        'EdgeColor', color_order(scaling_factor_index,:), ...
        'FaceColor', color_order(scaling_factor_index,:), ...
        'FaceAlpha', face_alphas(scaling_factor_index));
    xline(r_eff(scaling_factor_index), ...
        '-', ...
        'LineWidth', 2, ...
        'Color', color_order(scaling_factor_index,:));    
    box on;
    xlabel("$r [\mu$m$]$", 'Interpreter', 'latex');
    ylabel("probability", 'Interpreter', 'latex');
    pbaspect([5 4 1]);
    xticks(0:0.5:2.5);
    xlim([0 2.6]);
end

% save figure
print(gcf, '-dsvg', "scaled_axon_radius_distributions.svg");
