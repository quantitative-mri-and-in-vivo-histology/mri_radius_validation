function s = simulate_cylinder_signal_wpa(r, fiber_dir, bvec, bval, small_delta, big_delta, d_0)
% Simulates dMRI signal for cylinders using WPA.
%
%   Simulates the diffusion-weighted MRI signal for a set of impermeable 
%   cylinders with given radii, using the Weak Phase Approximation (WPA).
%
% USAGE:
%   s = simulate_cylinder_signal_wpa(r, fiber_dir, bvec, bval, small_delta, big_delta, d_0);
%
% INPUTS:
%   r           - (:,1) double  Cylinder radii.
%   fiber_dir   - (:,3) double  Fiber direction unit vectors.
%   bvec        - (:,3) double  Diffusion gradient directions.
%   bval        - (:,1) double  Diffusion b-values.
%   small_delta - (double)      Gradient pulse duration.
%   big_delta   - (double)      Gradient separation time.
%   d_0         - (double)      Free diffusivity.
%
% OUTPUT:
%   s - (array) Simulated diffusion signal using WPA.
%

arguments
    r (:,1) {mustBeNumeric,mustBeReal}
    fiber_dir (:,3) {mustBeNumeric,mustBeReal}
    bvec (:,3) {mustBeNumeric,mustBeReal}
    bval (:,1) {mustBeNumeric,mustBeReal}
    small_delta (:,1) {mustBeNumeric,mustBeReal}
    big_delta (:,1) {mustBeNumeric,mustBeReal}
    d_0 (1,1) {mustBeScalarOrEmpty}
end 

d_a_perp = 7 ./ 48 .* r'.^4 ./ (small_delta .* (big_delta - small_delta ./ 3) .* d_0);
d_a_perp = permute(shiftdim(d_a_perp, -1), [2 1 3]);
xi = bvec * fiber_dir';
s = exp(-bval .* d_a_perp .* (1 - xi.^2));

end
