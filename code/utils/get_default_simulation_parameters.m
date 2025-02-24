function params = get_default_simulation_parameters()
% Returns default parameters for simulations.
%
%   This function returns a struct containing default simulation settings 
%   such as noise realizations, random seed, generator type, radius 
%   approximation method, and bin edges.
%
% USAGE:
%   params = get_default_simulation_parameters();
%
% OUTPUT:
%   params - (struct) Default simulation parameters with the following fields:
%            * n_noise_realizations (int, default: 1000) - Number of noise realizations.
%            * random_seed (int, default: 0) - Random seed for reproducibility.
%            * random_generator (string, default: "Threefry") - Type of random generator.
%            * radius_approximation (string, default: "circular_equivalent") - 
%              Method for radius approximation.
%            * r_bin_edges (vector) - Bin edges for radius histograms.
%            * r_bin_centers (vector) - Corresponding bin centers.


    params = struct();
    params.n_noise_realizations = 1000;
    params.random_seed = 0;
    params.random_generator = "Threefry";
    params.radius_approximation = "circular_equivalent";
    params.r_bin_edges = unique([0:0.1:5 5:0.2:10 10:0.5:20])';
    params.r_bin_centers = params.r_bin_edges(1:end-1) ...
        + ((params.r_bin_edges(2:end) - params.r_bin_edges(1:end-1)) / 2);
end

