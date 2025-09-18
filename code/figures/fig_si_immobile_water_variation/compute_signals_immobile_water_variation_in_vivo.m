% ========================================================================
% Compute signals for discrete radii for clinical in-vivo dMRI protocols
% for various immobile water fraction values.
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

% simulation

% get in-vivo-like (scaled) axon radius distributions
[r_roi_counts, ...
    r_eff_reference] = get_histology_data(...
    fullfile(getenv("MRV_DATA_PATH")),...
    opt.sim_params.r_bin_edges,...
    opt.sim_params.radius_approximation,...
    "scaling_factor", opt.tissue_params.histology_scaling_factor);

tissue_params = opt.tissue_params;
tissue_params.f_im = 0.002;
sim_params = opt.sim_params;
snr_factors = opt.snr_factors;
r_eff_simulated = cell(1, length(opt.snr_factors));

states = cell(1, length(snr_factors));

f_im_step = 0.00005;
f_im_bin_centers = 0:f_im_step:tissue_params.f_im*2.5;
f_im_bin_edges = [f_im_bin_centers - f_im_step/2, f_im_bin_centers(end) + f_im_step/2];

signal_per_r = zeros(length(mr_protocols), ...
    length(f_im_bin_centers), ...
     sum(mr_protocol.num_gradient_directions_per_shell), ...
     length(sim_params.r_bin_centers));

for snr_baseline_index = 1:length(snr_factors)
    for f_im_index = 1:length(f_im_bin_centers)

  
        fprintf("Protocol: %3d/%3d; d_0: %04d/%04d\n", snr_baseline_index, length(snr_factors), f_im_index, length(f_im_bin_centers));

        tissue_params_adjusted = tissue_params;
        tissue_params_adjusted.f_im = f_im_bin_centers(f_im_index);
    
        % run simulations
        signal_per_r(snr_baseline_index,f_im_index,:,:) = simulate_signals_single_axon( ...
            sim_params.r_bin_centers, ...
            mr_protocols{snr_baseline_index}, ...
            tissue_params_adjusted, ...
            "signal_approximation", "matrix");
    end
end

r_bin_centers = sim_params.r_bin_centers;
r_bin_edges = sim_params.r_bin_edges;

save(fullfile(getenv("MRV_DATA_PATH"),...
    "simulations/signal_per_radius_in_vivo_clinical_var_immobile_water.mat"), ...
    "signal_per_r", ...
    "r_bin_centers", ...
    "r_bin_edges", ...
    "f_im_bin_centers", ...
    "f_im_bin_edges", ...
    "sim_params", ...
    "tissue_params", ...
    "mr_protocols");
