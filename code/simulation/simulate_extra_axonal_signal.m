function [s_extra, s_extra_par, s_extra_perp] = simulate_extra_axonal_signal(cylinder_dir, mr_protocol, d_e_par, d_e_perp)
% simulate_extra_axonal_signal Simulates the extra-axonal signal in diffusion MRI.
%
%   This function computes the extra-axonal signal contribution using a diffusion 
%   tensor model with parallel and perpendicular diffusivities.
%
% USAGE:
%   [s_extra, s_extra_par, s_extra_perp] = simulate_extra_axonal_signal( ...
%       cylinder_dir, mr_protocol, d_e_par, d_e_perp);
%
% INPUTS:
%   cylinder_dir - (Nx3 matrix) Fiber direction vectors.
%   mr_protocol  - (MrProtocol object) MRI protocol parameters.
%   d_e_par      - (scalar) Parallel diffusivity in the extra-axonal space.
%   d_e_perp     - (scalar) Perpendicular diffusivity in the extra-axonal space.
%
% OUTPUTS:
%   s_extra      - (MxN matrix) Total extra-axonal signal (M: gradient
%                   dirs)
%   s_extra_par  - (MxN matrix) Signal contribution from parallel diffusivity.
%   s_extra_perp - (MxN matrix) Signal contribution from perpendicular diffusivity.

    arguments
        cylinder_dir (:,3) {mustBeNumeric,mustBeReal}
        mr_protocol (1,1) MrProtocol
        d_e_par (1,1) {mustBeScalarOrEmpty}
        d_e_perp (1,1) {mustBeScalarOrEmpty}
    end

    xi_squared = (mr_protocol.bvec*cylinder_dir').^2;
    s_extra_perp = exp(-mr_protocol.bval_per_bvec.*d_e_perp.*(1-xi_squared));
    s_extra_par = exp(-mr_protocol.bval_per_bvec.*d_e_par.*xi_squared);
    s_extra = s_extra_par.*s_extra_perp;

end

