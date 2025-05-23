% ========================================================================
% Show effective radius correlation with histology for optimal 
% next-generation clinical scanner protocols
% ========================================================================

%% simulation
% load data
opt = load(fullfile(getenv("MRV_DATA_PATH"), ...
    "simulations/protocol_optimization_rois.mat"));
g_max_next_gen_clinical = 180;

% extract next-gen clinical scanner protocols for each SNR baseline level
mr_protocols = cell(1, length(opt.snr_factors));
signal_luts = cell(1, length(opt.snr_factors));
rand_substreams = zeros(1, length(opt.snr_factors));
for snr_baseline_index = 1:length(opt.snr_factors)
    
    % find protocol with optimal R
    roi_idx = abs(opt.results_table.g_max ...
        - g_max_next_gen_clinical) < 1e-1...
        & abs(opt.results_table.scaling_factor ...
        - opt.tissue_params.histology_scaling_factor) < 1e-3...
        & abs(opt.results_table.snr_factor ...
        - opt.snr_factors(snr_baseline_index)) < 1e-3;    
    
    roi_table = opt.results_table(roi_idx,:);
    [~, best_roi_idx] = max(roi_table.R, [], 'omitnan');
    
    % store protocol and signal lookup table
    best_protocol_index = roi_table.protocol_index(best_roi_idx);
    mr_protocol = opt.mr_protocols{best_protocol_index};
    mr_protocol.snr = roi_table.snr(best_roi_idx);
    signal_luts{snr_baseline_index} = ...
        squeeze(opt.signal_luts_intra_axonal(best_protocol_index, :, :));
    mr_protocols{snr_baseline_index} = mr_protocol;
end

%% simulation

% get in-vivo-like (scaled) axon radius distributions
[r_roi_counts, ...
    r_eff_reference] = get_histology_data(...
    fullfile(getenv("MRV_DATA_PATH")),...
    opt.sim_params.r_bin_edges,...
    opt.sim_params.radius_approximation,...
    "scaling_factor", opt.tissue_params.histology_scaling_factor);

tissue_params = opt.tissue_params;
sim_params = opt.sim_params;
snr_factors = opt.snr_factors;
r_eff_simulated = cell(1, length(opt.snr_factors));

states = cell(1, length(snr_factors));
for snr_baseline_index = 1:length(snr_factors)

    signals_per_roi = compute_ensemble_weighted_signal( ...
        sim_params.r_bin_centers, ...
        r_roi_counts, ...
        squeeze(signal_luts{snr_baseline_index}));

    rand_stream = RandStream(sim_params.random_generator, ...
        'Seed', sim_params.random_seed);
    r_eff_simulated{snr_baseline_index} = ...
    simulate_experimental_r_eff_estimation( ...
        signals_per_roi, ...
        mr_protocols{snr_baseline_index}, ...
        tissue_params.d_0, ...
        "noise_type", "gaussian", ...
        "noise_level", 1./mr_protocols{snr_baseline_index}.snr, ...
        "powder_average", "gaussian_ml", ...
        "f_im_estimate", tissue_params.f_im, ...
        "n_noise_realizations", sim_params.n_noise_realizations, ...
        "rand_stream", rand_stream);
end


%% plots
% set up figure
fig_handle  = figure;
set(fig_handle, get_default_figure_settings());
set(gcf,'units','centimeters','position',[0, 0, 18.4, 8.6]);
layout = tiledlayout(1, length(opt.snr_factors), ...
    'TileSpacing', 'tight', ...
    'Padding','tight');
plot_axes = gobjects(1, length(opt.snr_factors));
color_order = get_default_color_order();
snr_baseline_names = strings(1, length(opt.snr_factors));
for snr_factor_index = 1:length(opt.snr_factors)
    if opt.snr_factors(snr_factor_index) == 1
        snr_baseline_names(snr_factor_index) = "reference";
    else 
        snr_baseline_names(snr_factor_index) = ...
            sprintf("$%d \\%%$ increased", ...
            round(opt.snr_factors(snr_factor_index)*100-100));
    end
end

% plot simulated versus histological effective radii per SNR baseline level
for snr_baseline_index = 1:length(opt.snr_factors)

    % create inner layout to position title above legend
    inner_layout = tiledlayout(layout, 1, 1);
    inner_layout.Layout.Tile = snr_baseline_index;
    inner_layout.TileSpacing = "tight";
    inner_layout.Padding = "tight";
    title_text = "SNR baseline: " ...
        + snr_baseline_names(snr_baseline_index) + newline ...
        + sprintf("(SNR~$\\approx %d$)", ...
        int32(round(mr_protocols{snr_baseline_index}.snr)));
    title(inner_layout, title_text, ...
        "FontSize", 10, ...
        'interpreter', 'latex');   
    % plot correlations
    plot_axes(snr_baseline_index) = nexttile(inner_layout);
    rand_stream = RandStream(sim_params.random_generator, ...
        'Seed', sim_params.random_seed);
    plot_r_eff_correlation_ensemble(plot_axes(snr_baseline_index), ...
        r_eff_reference, ...
        r_eff_simulated{snr_baseline_index}, ...
        'color', color_order(snr_baseline_index,:), ...
        'xlim', [0 4], ...
        'ylim', [0 4], ...
        'rand_stream', rand_stream);
    xlabel('$r_{\mathrm{eff}}$ [$\mu$m] (histology)', ...
        'Interpreter', 'latex');
    ylabel(["$r_{\mathrm{eff}}$ [$\mu$m]", "(dMRI simulation)"], ...
        'Interpreter', 'latex');
    lgd = legend('Interpreter', 'latex', ...
        'Location', 'northoutside', ...
        'FontSize', 7);
    lgd.ItemTokenSize(1) = 18;
end

% save figure
print(gcf, '-dsvg', "optimized_correlation_next_gen_clinical.svg");
