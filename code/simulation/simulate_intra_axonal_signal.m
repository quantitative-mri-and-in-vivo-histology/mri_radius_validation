function s = simulate_intra_axonal_signal(r, fiber_dir, mr_protocol, d_0, d_a_par, model)
% Simulates intra-axonal dMRI signal using various models.
%
%   Computes the intra-axonal diffusion-weighted MRI signal using different 
%   modeling approaches: Gaussian Phase Approximation (GPA), Weak Phase 
%   Approximation (WPA), or the Matrix Method.
%
% USAGE:
%   s = simulate_intra_axonal_signal(r, fiber_dir, mr_protocol, d_0, d_a_par, model);
%
% INPUTS:
%   r           - (:,1) double    Cylinder radii.
%   fiber_dir   - (:,3) double    Fiber direction unit vectors.
%   mr_protocol - (1,1) struct    MR protocol containing bvec, bval, and timing.
%   d_0         - (double)        Free diffusivity.
%   d_a_par     - (double)        Axial diffusivity along the fiber.
%   model       - (string)        Model type: "gpa", "wpa", or "matrix".
%
% OUTPUT:
%   s - (array) Simulated intra-axonal diffusion signal.
%

arguments
    r (:,1) {mustBeNumeric,mustBeReal}
    fiber_dir (:,3) {mustBeNumeric,mustBeReal}
    mr_protocol (1,1) MrProtocol
    d_0 (1,1) {mustBeScalarOrEmpty}
    d_a_par (1,1) {mustBeScalarOrEmpty}
    model (1,1) string
end 

if model == "gpa"
    s = simulate_intra_axonal_signal_gpa(r, fiber_dir, ...
        mr_protocol.bvec, mr_protocol.bval_per_bvec, ...
        mr_protocol.gradient_duration_per_bvec, ...
        mr_protocol.gradient_separation_per_bvec, d_0, d_a_par);
    
elseif model == "wpa"
    s = simulate_intra_axonal_signal_wpa(r, fiber_dir, ...
        mr_protocol.bvec, mr_protocol.bval_per_bvec, ...
        mr_protocol.gradient_duration_per_bvec, ...
        mr_protocol.gradient_separation_per_bvec, d_0, d_a_par);
    
elseif model == "matrix"
    s = simulate_intra_axonal_signal_matrix_method(r, fiber_dir, ...
        mr_protocol.bvec, mr_protocol.bval_per_bvec, ...
        mr_protocol.gradient_duration_per_bvec, ...
        mr_protocol.gradient_separation_per_bvec, d_0, d_a_par);
    
else
    error("Unknown model. Choose 'gpa', 'wpa', or 'matrix'.");
end

end
