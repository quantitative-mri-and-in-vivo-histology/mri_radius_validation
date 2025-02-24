function r_eff = compute_r_eff_from_distribution(varargin)
% Computes r_eff from radii or histogram data.
%
%   Computes the effective radius (r_eff) from either a direct set of radii 
%   or from histogram data (bin centers and bin counts).
%
% USAGE:
%   r_eff = compute_r_eff_from_distribution(radii)
%   r_eff = compute_r_eff_from_distribution('BinCenters', bin_centers, ...
%           'BinCounts', bin_counts)
%
% INPUTS:
%   radii        - (vector, optional) A list of radii values.
%   'BinCenters' - (vector, optional) Centers of histogram bins.
%   'BinCounts'  - (vector, optional) Corresponding counts per bin.
%
%   * Either provide a direct list of radii OR specify both 'BinCenters' 
%   and 'BinCounts'.
%
% OUTPUT:
%   r_eff - (double) Effective radius computed as:
%           r_eff = (sum(R.^6) / sum(R.^2))^(1/4)

    if nargin == 1 && isnumeric(varargin{1}) % Case 1: Direct radii input
        radii = varargin{1};
        r_eff = (sum(radii.^6) ./ sum(radii.^2)).^(1/4);
    
    elseif nargin == 4 && ischar(varargin{1}) && strcmp(varargin{1}, ...
                'BinCenters') ...
                      && ischar(varargin{3}) && strcmp(varargin{3}, ...
                 'BinCounts')
        % Case 2: Histogram data input
        bin_centers = varargin{2};
        bin_counts = varargin{4};
        r_eff = (sum(bin_centers.^6 .* bin_counts) ./ sum(bin_centers.^2 .* bin_counts)).^(1/4);

    else
        error(['Invalid input. Use either:\n' ...
               '  compute_r_eff_from_distribution(radii)\n' ...
               '  compute_r_eff_from_distribution(''BinCenters'', bin_centers, ''BinCounts'', bin_counts).']);
    end
end
