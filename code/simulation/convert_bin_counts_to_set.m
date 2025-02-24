function values = convert_bin_counts_to_set(counts, bin_centers)
% convert_bin_counts_to_set Expands histogram bin counts into a value set.
%
%   Converts histogram bin counts and their corresponding bin centers into 
%   a set of individual values.
%
% USAGE:
%   values = convert_bin_counts_to_set(counts, bin_centers);
%
% INPUTS:
%   counts      - (vector) Number of occurrences in each histogram bin.
%   bin_centers - (vector) Center values of the histogram bins.
%
% OUTPUT:
%   values - (vector) Expanded dataset where each bin center is repeated 
%            according to its count.
%

values = zeros(sum(counts(:)), 1);
n_values = 1;
for i = 1:length(bin_centers)
    if counts(i) > 0                  
        values(n_values:n_values + counts(i) - 1) = bin_centers(i);
        n_values = n_values + counts(i);
    end
end

end
