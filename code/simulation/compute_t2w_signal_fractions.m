function [f_a, f_e] = compute_t2w_signal_fractions(t_e, f_0, t_2_a, t_2_e, f_im)
% Computes signal fractions for a T2-weighted signal.
%
%   Computes the signal fractions of axonal (f_a) and extracellular (f_e) 
%   compartments in a multi-compartment T2-weighted model, accounting for an 
%   immobile compartment (f_im).
%
% USAGE:
%   [f_a, f_e] = compute_t2w_signal_fractions(t_e, f_0, t_2_a, t_2_e, f_im);
%
% INPUTS:
%   t_e   - (double) Echo time (TE) in seconds.
%   f_0   - (double) Initial fraction of the axonal compartment.
%   t_2_a - (double) T2 relaxation time of the axonal compartment.
%   t_2_e - (double) T2 relaxation time of the extracellular compartment.
%   f_im  - (double) Fraction of the immobile compartment.
%
% OUTPUTS:
%   f_a - (double) Signal fraction from the axonal compartment.
%   f_e - (double) Signal fraction from the extracellular compartment.

% Total signal denominator
s_total = f_0 * exp(-t_e / t_2_a) + (1 - f_0 - f_im) .* exp(-t_e / t_2_e) + f_im;

% Axonal signal fraction
f_a = (f_0 * exp(-t_e / t_2_a)) ./ s_total;

% Extracellular signal fraction
f_e = ((1 - f_0 - f_im) .* exp(-t_e / t_2_e)) ./ s_total;
end
