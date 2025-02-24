function s = compute_ensemble_weighted_signal(r_bin_centers, r_counts, s_per_r_bin)
    
% Computes ensemble-weighted signal based on signal lookup table.
%
% INPUTS:
%   r_bin_centers   - (Nx1 matrix) Radii bins.
%   r_counts        - (MxN matrix) Counts per radii bin. (M = samples)
%   s_per_bin       - (LxN matrix) Signals per radii bin (L = grad. dirs)
%
% OUTPUT:
%   s - (LxM matrix) Ensemble-weighted signals.

    arguments
        r_bin_centers (:,1) {mustBeNumeric,mustBeReal}
        r_counts (:,:) {mustBeNumeric,mustBeReal}
        s_per_r_bin (:,:) {mustBeNumeric,mustBeReal}
    end 
    
    assert(size(r_bin_centers, 1) == size(r_counts, 2));
    assert(size(r_bin_centers, 1) == size(s_per_r_bin, 2))

    norm_term = r_counts.*r_bin_centers'.^2;
    s = tensorprod(norm_term, s_per_r_bin, 2, 2)./sum(norm_term,2);
   
    if ndims(s) > 2
        perm_order = [2 1 3:1:ndims(s)];
    else
        perm_order = [2 1];
    end

    s = permute(s, perm_order);
end
