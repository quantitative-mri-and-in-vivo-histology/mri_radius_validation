function s = simulate_intra_axonal_signal_matrix_method(r, fiber_dir, bvec, bval, small_delta, big_delta, d_0, d_a_par)
% Simulates intra-axonal dMRI signal using the Matrix Method.
%
%   Computes the intra-axonal diffusion-weighted MRI signal using the Matrix 
%   Method, accounting for both parallel and perpendicular diffusion components. 
%   The perpendicular signal is simulated using a numerical approach with waveform 
%   discretization and time evolution matrices.
%
% USAGE:
%   s = simulate_intra_axonal_signal_matrix_method(r, fiber_dir, bvec, bval, ...
%       small_delta, big_delta, d_0, d_a_par);
%
% INPUTS:
%   r           - (:,1) double  Cylinder radii (in µm).
%   fiber_dir   - (:,3) double  Fiber direction unit vectors.
%   bvec        - (:,3) double  Diffusion gradient directions.
%   bval        - (:,1) double  Diffusion b-values (in s/mm²).
%   small_delta - (double)      Gradient pulse duration (in ms).
%   big_delta   - (double)      Gradient separation time (in ms).
%   d_0         - (double)      Free diffusivity (in µm²/ms).
%   d_a_par     - (double)      Axial diffusivity along the fiber (in µm²/ms).
%
% OUTPUT:
%   s - (array) Simulated intra-axonal diffusion signal.
%

arguments
    r (:,1) {mustBeNumeric,mustBeReal}
    fiber_dir (:,3) {mustBeNumeric,mustBeReal}
    bvec (:,3) {mustBeNumeric,mustBeReal}
    bval (:,1) {mustBeNumeric,mustBeReal}
    small_delta (:,1) {mustBeNumeric,mustBeReal}
    big_delta (:,1) {mustBeNumeric,mustBeReal}
    d_0 (1,1) {mustBeScalarOrEmpty}
    d_a_par (1,1) {mustBeScalarOrEmpty}
end 

tau = 1E-5; % Time interval for waveform discretization.

% Convert units to SI
small_delta = small_delta .* 1e-3; % ms → s
big_delta = big_delta .* 1e-3;     % ms → s
d_0 = d_0 * 1e-9;                 % µm²/ms → m²/s
d_a_par = d_a_par * 1e-9;         % µm²/ms → m²/s
r = r * 1e-6;                     % µm → m
bval = bval * 1e9;                % s/mm² → s/m²

gyroMagnRatio = 2.67513 * 10^(8); % Gyromagnetic ratio (rad/T·s)
g = sqrt(bval ./ (small_delta.^2 .* (gyroMagnRatio).^2 .* (big_delta - small_delta / 3)));

s_perp = simulate_cylinder_signal_matrix_method(r, ...
    fiber_dir, g, bvec, unique(small_delta), unique(big_delta), d_0, tau);

s_perp = permute(s_perp, [3 2 1]);
xi = (bvec * fiber_dir');
s_par = exp(-bval .* d_a_par .* xi.^2);
s = abs(s_par .* s_perp);

end
