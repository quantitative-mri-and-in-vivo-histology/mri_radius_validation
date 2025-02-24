function [r_bin_counts, r_eff] = get_histology_data(data_path, r_bin_edges, radius_approximation, options)
% Loads histological axon radius distributions and computes effective radii.
%
%   This function loads axon radius histograms from histology data and 
%   computes the effective radius (r_eff) based on the chosen radius 
%   approximation method.
%
% USAGE:
%   [r_bin_counts, r_eff] = get_histology_data(data_path, r_bin_edges, ...
%       "circular_equivalent", "scaling_factor", 1.2);
%
% INPUTS:
%   data_path           - (string) Path to the histology raw data directory.
%   r_bin_edges         - ((N+1)x1 vector) Bin edges for radius histogram (N: # bins).
%   radius_approximation - (string) Method for approximating radius:
%                          * "minor_axis" for minor axis length.
%                          * "circular_equivalent" for circular equivalent radius.
%   options.scaling_factor - (scalar, default: 1.0) Factor to scale radii.
%
% OUTPUTS:
%   r_bin_counts - (MxN matrix) Histogram counts of axon radii per
%                   histology section. (M: samples)
%   r_eff        - (Mx vector) Effective radius computed from distributions.

    arguments
        data_path (:,1) string
        r_bin_edges (:,1) {mustBeReal, mustBeNumeric}
        radius_approximation (1,1) string
        options.scaling_factor (1,1) {mustBeReal, mustBeNumeric} = 1.0
    end 

    histology_data_dir = fullfile(data_path, "histology", "rawdata");
    r_bin_centers_org = readtable(fullfile( ...
        histology_data_dir, "desc-binCenters_radii.tsv"), ...
        "FileType", "text", ...
        "Delimiter", "\t");
    r_bin_centers_org = r_bin_centers_org.Variables;

    r_bin_centers = r_bin_edges(1:end-1) ...
        + ((r_bin_edges(2:end) - r_bin_edges(1:end-1)) / 2);

    if radius_approximation == "minor_axis"
        r_bin_counts_org = readtable(fullfile( ...
            histology_data_dir, "desc-countsMinorAxis_radii.tsv"), ...
            "FileType", "text", ...
            "Delimiter", "\t");
        r_bin_counts_org = r_bin_counts_org.Variables;
    elseif radius_approximation == "circular_equivalent"
        r_bin_counts_org = readtable(fullfile( ...
            histology_data_dir, "desc-countsCircularEq_radii.tsv"), ...
            "FileType", "text", ...
            "Delimiter", "\t");
        r_bin_counts_org = r_bin_counts_org.Variables;
    else
        error("Unknwon radius type.")
    end

    r_bin_counts = zeros(size(r_bin_counts_org,1), length(r_bin_centers));
    r_eff = zeros(1, size(r_bin_counts_org,1));

    for section_index = 1:size(r_bin_counts_org,1); 
        axon_radii = convert_bin_counts_to_set( ...
            r_bin_counts_org(section_index,:), ...
            r_bin_centers_org);
        axon_radii_scaled = axon_radii*options.scaling_factor;
        r_section_counts_scaled = histcounts( ...
            axon_radii_scaled, ...
            r_bin_edges);
        r_bin_counts(section_index, :) = r_section_counts_scaled;
        r_eff(section_index) = compute_r_eff_from_distribution( ...
            'BinCenters', r_bin_centers', ...
            'BinCounts', r_section_counts_scaled);
    end

end