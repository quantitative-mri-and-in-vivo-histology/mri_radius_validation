function s = simulate_cylinder_signal_gpa(r, fiber_dir, bvec, bval, delta, Delta, d_0)
% Simulates dMRI signal for cylinders using GPA.
%
%   Simulates the diffusion-weighted MRI signal for a set of impermeable 
%   cylinders with given radii, using the Gaussian Phase Approximation (GPA).
%
% USAGE:
%   s = simulate_cylinder_signal_gpa(r, fiber_dir, bvec, bval, delta, Delta, d_0);
%
% INPUTS:
%   r         - (:,1) double  Cylinder radii.
%   fiber_dir - (:,3) double  Fiber direction unit vectors.
%   bvec      - (:,3) double  Diffusion gradient directions.
%   bval      - (:,1) double  Diffusion b-values.
%   delta     - (:,1) double  Gradient pulse durations.
%   Delta     - (:,1) double  Gradient separation times.
%   d_0       - (1,1) double  Free diffusivity.
%
% OUTPUT:
%   s - (array) Simulated diffusion signal.
%

arguments
    r (:,1) {mustBeNumeric,mustBeReal}
    fiber_dir (:,3) {mustBeNumeric,mustBeReal}
    bvec (:,3) {mustBeNumeric,mustBeReal}
    bval (:,1) {mustBeNumeric,mustBeReal}
    delta (:,1) {mustBeNumeric,mustBeReal}
    Delta (:,1) {mustBeNumeric,mustBeReal}
    d_0 (1,1) {mustBeScalarOrEmpty}
end 

td = r'.^2 / d_0;
bardelta = delta ./ td;
barDelta = Delta ./ td;

N = 10; 
alpha = [1.8412 5.3314 8.5363 11.7060 14.8636 18.0155 21.1644 ...
         24.3113 27.4571 30.6019 33.7462 36.8900 40.0334 43.1766 ...
         46.3196 49.4624 52.6050 55.7476 58.8900 62.0323];
alpha = alpha(1:N);
alpha = shiftdim(alpha, -1);

s = (2 ./ (alpha.^6 .* (alpha.^2 - 1))) .* (-2 + 2 .* alpha.^2 .* bardelta + ...
     2 .* (exp(-alpha.^2 .* bardelta) + exp(-alpha.^2 .* barDelta)) - ...
     exp(-alpha.^2 .* (bardelta + barDelta)) - exp(-alpha.^2 .* (barDelta - bardelta)));

s = sum(s, 3);
q = 1 ./ delta .* sqrt(bval ./ (Delta - delta ./ 3));
s = -s .* d_0 .* q.^2 .* td.^3;
s = permute(shiftdim(s, -1), [2 1 3]);
s = exp(s .* (1 - (bvec * fiber_dir').^2));

end
