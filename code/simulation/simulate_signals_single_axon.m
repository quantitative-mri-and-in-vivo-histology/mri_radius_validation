function [s_total, s_intra, s_extra] = simulate_signals_single_axon(...
        r, mr_protocol, tissue_params, options)
% Simulates dMRI signals for a single axon.
%
%   Computes intra-axonal, extra-axonal, and total diffusion-weighted MRI 
%   signals for a single axon, incorporating T2-weighted signal fractions 
%   and integrating over fiber orientation dispersion using Lebedev quadrature.
%
% USAGE:
%   [s_total, s_intra, s_extra] = simulate_signals_single_axon(...
%       r, mr_protocol, tissue_params, options);
%
% INPUTS:
%   r           - (N,1) double  Cylinder radii.
%   mr_protocol - (1,1) struct  MR protocol containing acquisition parameters.
%   tissue_params - (1,1) struct  Tissue parameters affecting signal fractions.
%   options.signal_approximation - (string) Signal computation method ("matrix", default).
%   options.n_lebedev_directions - (scalar) Number of Lebedev quadrature directions (default: 590).
%   options.intra_axonal_signal_table - (array) Precomputed intra-axonal signal lookup table.
%   options.cylinder_dir - (1,3) double  Preferred cylinder direction (default: [0 0 1]).
%
% OUTPUTS:
%   s_total - (MxN matrix) Total diffusion-weighted MRI signal (M: gradient
%               dirs)
%   s_intra - (MxN matrix) Intra-axonal diffusion signal.
%   s_extra - (Mx1 matrix) Extra-axonal diffusion signal.
%

arguments
    r (:,1) {mustBeNumeric,mustBeReal}
    mr_protocol (1,1) MrProtocol
    tissue_params (1,1) struct
    options.signal_approximation (1,1) string = "matrix";
    options.n_lebedev_directions (1,1) {mustBeScalarOrEmpty} = 590
    options.intra_axonal_signal_table {mustBeNumeric,mustBeReal} = []
    options.cylinder_dir (1,3) {mustBeNumeric,mustBeReal} = [0 0 1];
end 

[f_a, f_e] = compute_t2w_signal_fractions(mr_protocol.t_e, ...
    tissue_params.f_0, ...
    tissue_params.t_2_a, ...
    tissue_params.t_2_e, ...
    tissue_params.f_im);

% Compute fiber directions and weights for Lebedev integration
lebedev_sphere = compute_lebedev_sphere(options.n_lebedev_directions);
lebedev_dirs = [lebedev_sphere.x lebedev_sphere.y lebedev_sphere.z];
lebedev_weights = lebedev_sphere.w;
watson_factor = (4 * pi * hypergeom((1/2), (3/2), tissue_params.neurite_dispersion))^-1;
watson_weights = watson_factor * exp(tissue_params.neurite_dispersion * (options.cylinder_dir * lebedev_dirs')'.^2);
combined_weights = lebedev_weights .* watson_weights;

% Compute extra-axonal signal
if (f_e > 0 && ~isempty(tissue_params.d_e_par) && isfinite(tissue_params.d_e_par) ...
        && ~isempty(tissue_params.d_e_perp) && isfinite(tissue_params.d_e_perp))
    s_extra_per_fiber_dir = simulate_extra_axonal_signal(lebedev_dirs, ...
        mr_protocol, ...
        tissue_params.d_e_par, ...
        tissue_params.d_e_perp);
else
    s_extra_per_fiber_dir = zeros(size(mr_protocol.bvec, 1), size(lebedev_dirs, 1));
end
s_extra = squeeze(f_e .* tensorprod(s_extra_per_fiber_dir, combined_weights, 2, 1));

% Compute intra-axonal signal
if (f_a > 0 && ~isempty(tissue_params.d_a_par) && isfinite(tissue_params.d_a_par))
    s_intra_per_fiber_dir = simulate_intra_axonal_signal(r, lebedev_dirs, ...
        mr_protocol, tissue_params.d_0, tissue_params.d_a_par, options.signal_approximation);
end
s_intra = f_a .* tensorprod(s_intra_per_fiber_dir, combined_weights, 2, 1);
s_intra = squeeze(permute(s_intra, [1 3 2]));

% Combine signals
s_total = s_intra + s_extra + tissue_params.f_im;

end
