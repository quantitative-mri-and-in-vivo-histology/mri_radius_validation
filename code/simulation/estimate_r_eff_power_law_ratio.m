function [r_eff, D_perp, log_arg] = estimate_r_eff_power_law_ratio(s, b, delta, Delta, D0)
% Estimates effective radius using a power-law ratio of powder-averaged
% signals over two diffusion shells. The approach was described by 
% Pizzolato et. al (https://doi.org/10.1016/j.media.2023.102767)
%
% USAGE:
%   [r_eff, D_perp, log_arg] = ...
% estimate_r_eff_power_law_ratio(s, b, delta, Delta, D0);
%
% INPUTS:
%   s      - (array) Powder-averaged signals for different b-values.
%   b      - (array) Diffusion b-values.
%   delta  - (scalar) Gradient pulse duration.
%   Delta  - (scalar) Gradient separation time.
%   D0     - (scalar) Free diffusivity.
%
% OUTPUTS:
%   r_eff   - (array) Estimated effective axon radius.
%   D_perp  - (array) Estimated perpendicular diffusivity.
%   log_arg - (array) Argument used in the logarithmic term for D_perp computation.
%

s_ratio = s(1,:) ./ s(end,:);
s_ratio_invalid_idx = isinf(s_ratio) | isnan(s_ratio) | s_ratio <= 0;
s_ratio(s_ratio_invalid_idx) = 1;

log_arg = s_ratio .* sqrt(b(1) ./ b(end));
D_perp = log(log_arg) ./ (b(end) - b(1));

r_eff = nan(size(D_perp));
r_eff(D_perp > 0) = (48 / 7 .* delta .* (Delta - delta / 3) ...
    .* D0 .* D_perp(D_perp > 0)).^(1/4);

end
