function s = simulate_intra_axonal_signal_wpa(r, fiber_dir, bvec, bval, small_delta, big_delta, d_0, d_a_par)
% Simulates intra-axonal dMRI signal using WPA.
%
%   Computes the intra-axonal diffusion-weighted MRI signal using the 
%   Weak Phase Approximation (WPA). The total signal is modeled as the 
%   product of perpendicular and parallel diffusion components.
%
% USAGE:
%   s = simulate_intra_axonal_signal_wpa(r, fiber_dir, bvec, bval, ...
%       small_delta, big_delta, d_0, d_a_par);
%
% INPUTS:
%   r           - (:,1) double  Cylinder radii.
%   fiber_dir   - (:,3) double  Fiber direction unit vectors.
%   bvec        - (:,3) double  Diffusion gradient directions.
%   bval        - (:,1) double  Diffusion b-values.
%   small_delta - (double)      Gradient pulse duration.
%   big_delta   - (double)      Gradient separation time.
%   d_0         - (double)      Free diffusivity.
%   d_a_par     - (double)      Axial diffusivity along the fiber.
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

assert(size(bval,1) == size(bvec,1));

xi = (bvec * fiber_dir');
s_perp = simulate_cylinder_signal_wpa(r, fiber_dir, bvec, bval, small_delta, big_delta, d_0);
s_par = exp(-bval .* d_a_par .* xi.^2);
s = s_par .* s_perp;

end
