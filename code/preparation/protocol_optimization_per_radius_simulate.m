% ========================================================================
% Compute signal LUT for various in-vivo protocol candidates
% ========================================================================

% load in-vivo parameters
sim_params = get_default_simulation_parameters();
tissue_params = jsondecode(fileread(...
    "../../parameters/tissue_params_in_vivo.json"));
mr_protocol_reference = MrProtocol.fromJsonFile(....
    "../../parameters/mr_protocol_in_vivo_experimental.json");

% define physics constants
gyroMagnRatio =  267.513*10^(-6);

% define fixed parameters
b_min = 6;
t_tr = mr_protocol_reference.t_e ...
    - mr_protocol_reference.gradient_duration_per_shell(1) ...
    - mr_protocol_reference.gradient_separation_per_shell(1);

% compute maximum gradient amplitude of existing scanners (assuming 90
% percent of nominal values)
g_nominal_factor = 0.9;
g_max_scanners_theoretical = [80 200 300 500]';
g_max_scanners = g_nominal_factor.*g_max_scanners_theoretical;

% define search grid
b_max = 800;
t_e_max = 400;
g_max = linspace(min(10), max(600), 15)';
g_max = sort([g_max; g_nominal_factor .* g_max_scanners_theoretical]);
small_delta = unique([2:2:40 40:4:60])';
big_delta = unique([2:2:60 60:4:80]);

[g_high_grid, small_delta_grid, big_delta_grid] = ndgrid( ...
    g_max, small_delta, big_delta);

g_min_grid = nan(size(g_high_grid));
b_max_grid = nan(size(g_high_grid));
t_e_grid = nan(size(g_high_grid));

% compute dependent variables
for i = 1:numel(g_high_grid)
    t_e_grid(i) = big_delta_grid(i) + small_delta_grid(i) + t_tr;
    g_min_grid(i) = sqrt(b_min./(small_delta_grid(i).^2 ...
        .*gyroMagnRatio.^2.*(big_delta_grid(i)-small_delta_grid(i)/3)));
    b_max_grid(i) = (g_high_grid(i).*gyroMagnRatio ...
        .*small_delta_grid(i)).^2.*(big_delta_grid(i) ...
        -small_delta_grid(i)/3);
end

% filter out implausible protocols
nan_idx = (imag(g_min_grid) > 0) ...
    | (imag(b_max_grid) > 0) ...
    | (b_max_grid < b_min) ...
    | (b_max_grid > b_max) ...
    | (t_e_grid > t_e_max) ...
    | (small_delta_grid + 4 >= big_delta_grid) ...
    | (g_high_grid <= g_min_grid) ...
    | isnan(g_min_grid) ...
    | isnan(b_max_grid) ...
    | isnan(t_e_grid);
valid_idx = find(~nan_idx);

% create protocols from parameters
mr_protocols = cell(length(valid_idx), 1);
for good_config_index = 1:length(valid_idx)
    grid_index = valid_idx(good_config_index);
    
    mr_protocol = MrProtocol( ...
        [g_min_grid(grid_index); g_high_grid(grid_index)], ...
        repmat(small_delta_grid(grid_index), ...
        mr_protocol_reference.num_shells, 1), ...
        repmat(big_delta_grid(grid_index), ...
        mr_protocol_reference.num_shells, 1), ...
        [b_min; b_max_grid(grid_index)], ...
        mr_protocol_reference.bvec_per_shell, ...
        t_e_grid(grid_index));
    snr_scaling_factor = ...
        (tissue_params.f_0*(exp(-mr_protocol.t_e/tissue_params.t_2_a)) ...
        + (1-tissue_params.f_0)...
        *(exp(-mr_protocol.t_e/tissue_params.t_2_e)))...
        /(tissue_params.f_0...
        *exp(-mr_protocol_reference.t_e/tissue_params.t_2_a) ...
        + (1-tissue_params.f_0) ... 
        .*exp(-mr_protocol_reference.t_e/tissue_params.t_2_e));
    mr_protocol.snr = snr_scaling_factor.*mr_protocol_reference.snr;
    mr_protocols{good_config_index} = mr_protocol;
end

% compute signal LUT per protocol
signal_luts_intra_axonal = zeros(length(mr_protocols), ...
    length(mr_protocol_reference.bval_per_bvec), ...
    length(sim_params.r_bin_centers));
data_queue = parallel.pool.DataQueue;
wbar = waitbar(0, 'Computing signal LUTs...');
n = length(mr_protocols);
parallel_waitbar_update(n, wbar);
afterEach(data_queue, @parallel_waitbar_update);
parfor valid_index = 1:n
    signal_per_r = simulate_signals_single_axon( ...
        sim_params.r_bin_centers, ...
        mr_protocols{valid_index}, ...
        tissue_params, ...
        "signal_approximation", "matrix");
    signal_luts_intra_axonal(valid_index,:,:) = signal_per_r;
    send(data_queue, 1);
end
close(wbar);

% save protocols and corresponding signal LUTs
save(fullfile(getenv("MRV_DATA_PATH"), ...
    "simulations/protocol_optimization_per_radius.mat"), ...
    'mr_protocols', ...
    'signal_luts_intra_axonal', ...
    'mr_protocol_reference', ...
    'g_max_scanners', ...
    'tissue_params', ...
    'sim_params');
