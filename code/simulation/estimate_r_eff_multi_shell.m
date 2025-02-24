function [r_eff, beta] = estimate_r_eff_multi_shell(s, b, delta, Delta, D0)    
% estimate_r_eff_multi_shell Estimates r_eff from multi-shell dMRI data.
%
%   Estimates the effective axon radius (r_eff), beta parameter
%   using diffusion MRI data acquired with multiple b-values.
%
% USAGE:
%   [r_eff, beta, D_perp] = estimate_r_eff_multi_shell(s, b, delta, Delta, D0);
%
% INPUTS:
%   s      - (NxM double) Signal intensities (N samples, M b-values).
%   b      - (1xM double) Diffusion b-values.
%   delta  - (1xM double) Gradient pulse durations.
%   Delta  - (1xM double) Gradient separation times.
%   D0     - (double) Diffusivity of the axoplasm.
%
% OUTPUTS:
%   r_eff  - (Nx1 double) Estimated effective axon radii.
%   beta   - (Nx1 double) Estimated beta. (slope of the signal scaling as 
%               a function of 1/sqrt. ~ f/sqrt(D_{a,perp}))
% 
% values from fitting.
%

g = sqrt(b ./ (delta.^2 .* (267.513 * 10^(-6)).^2 .* (Delta - delta/3)));
r_eff = nan(size(s,1), 1);
beta = nan(size(s,1), 1);

parfor i = 1:size(s,1)
    values = s(i,:);
    if all(values)
        [r_eff(i), beta(i)] = estimate_r_eff_non_linear_fit( ...
            delta, Delta, g, values', "Neuman", D0);
    end
end

invalid_idx = r_eff < 0.1;
r_eff(invalid_idx) = nan;
beta(invalid_idx) = nan;

end
