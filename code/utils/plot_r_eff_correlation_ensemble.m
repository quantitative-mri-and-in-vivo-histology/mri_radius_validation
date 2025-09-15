function ax_handle = plot_r_eff_correlation_ensemble(ax_handle, r_eff_reference, r_eff_estimated, options)
% Plots correlation between histological and estimated r_eff.
%
%   This function visualizes the correlation between histological reference 
%   effective axon radii (r_eff) and their corresponding estimated values 
%   from dMRI data, optionally including statistical metrics, regression lines, 
%   and confidence intervals.
%
% USAGE:
%   ax_handle = plot_r_eff_correlation_ensemble(ax_handle, ...
%       r_eff_reference, r_eff_estimated, options);
%
% INPUTS:
%   ax_handle       - (axes handle) Target axis for plotting.
%   r_eff_reference - (Nx1 vector) Reference r_eff values from histology.
%   r_eff_estimated - (NxL matrix) Estimated r_eff values, where L is the 
%                     number of noise realizations.
%   options         - (struct) Plot customization options:
%                     * marker (string, default: "o") - Marker style.
%                     * marker_size (int, default: 3) - Marker size.
%                     * regression_line_style (string, default: ":") - Line style for regression.
%                     * plot_regression (bool, default: true) - Whether to plot regression.
%                     * plot_error_metrics (bool, default: true) - Whether to annotate correlation stats.
%                     * plot_unity_line (bool, default: true) - Whether to plot unity line.
%                     * plot_confidence_interval (bool, default: true) - Whether to show 95% CI.
%                     * color (1x3 vector, default: pink RGB) - Plot color.
%                     * legend_prefix (string, default: "") - Prefix for legend entries.
%                     * xlim, ylim (1x2 vectors, default: auto) - Axes limits.
%                     * rand_stream (RandStream, default: global) - Random number generator.
%
% OUTPUT:
%   ax_handle - (axes handle) Handle to the modified plot.
%

arguments
    ax_handle (1,1)
    r_eff_reference (:,:) {mustBeNumeric}
    r_eff_estimated (:,:,:) {mustBeNumeric}
    options.marker (1,1) string = "o"
    options.marker_size (1,1) int32 = 3
    options.regression_line_style (1,1) string = ":"
    options.plot_regression (1,1) logical = true
    options.plot_error_metrics (1,1) logical = true
    options.plot_unity_line (1,1) logical = true
    options.plot_confidence_interval (1,1) logical = true
    options.color (1,3) {mustBeNumeric} = 1./255.*[255, 102, 204]
    options.legend_prefix (1,1) string = ""
    options.xlim (1,2) {mustBeNumeric} = []
    options.ylim (1,2) {mustBeNumeric} = []
    options.rand_stream = []
end

if isempty(options.rand_stream)
    rand_stream = RandStream.getGlobalStream();
else
    rand_stream = options.rand_stream;
end

assert(ndims(r_eff_estimated) >= 2 & ndims(r_eff_estimated) <= 3); 
ci_face_alpha = 0.2;
has_noise_realizations = sum(size(r_eff_estimated) > 1) > 1;

axes(ax_handle);     
hold on;   

% Determine axis limits
if ~isempty(options.xlim) && ~isempty(options.ylim)
    xlims = options.xlim;
    ylims = options.ylim;
else
    r_eff_max = ceil(max([r_eff_reference(:); r_eff_estimated(:)]));
    xlims = [0 r_eff_max];
    ylims = xlims;
end
xtick_vals = xlims(1):1:xlims(2);
ytick_vals = ylims(1):1:ylims(2);

% Compute correlation statistics
[r_eff_median_values, r_eff_conf_values, correlation_r, correlation_p, ...
 linear_model, nrmse, fit_success_ratio] = compute_r_eff_correlation_metrics( ...
    squeeze(r_eff_reference)', squeeze(r_eff_estimated), 'rand_stream', rand_stream);

% Plot confidence intervals
if options.plot_confidence_interval
    [~, sorted_idx] = sort(r_eff_reference);
    fill([r_eff_reference(sorted_idx) fliplr(r_eff_reference(sorted_idx))], ...
        [r_eff_conf_values(sorted_idx,1)' flipud(r_eff_conf_values(sorted_idx,2))'], ...
        options.color, 'FaceAlpha', ci_face_alpha, 'EdgeColor', [1 1 1], 'HandleVisibility', 'off');
end

% Prepare legend string
if has_noise_realizations 
    legend_str = string(options.legend_prefix) + "ROI medians + 95$\%$ CI";
else
    legend_str = string(options.legend_prefix) + "ROIs";
end

if options.plot_error_metrics    
    if correlation_p < 1e-3
        p_str = sprintf("p < 1e^{-3}");
    else
        p_str = sprintf("p = %.1f e^{%d}",  ...
            correlation_p / 10^floor(log10(correlation_p)), floor(log10(correlation_p)));
    end

    if has_noise_realizations
        legend_str = legend_str + newline ...
            + sprintf("$\\quad(R = %.2f, %s$" + newline ...
            + "$\\quad \\mathrm{NRMSE} = %.1f \\%%, S = %.1f \\%%$)", ...
            correlation_r, p_str, 100*nrmse, 100*fit_success_ratio);
    else
        if (string(options.legend_prefix) ~= "")
            legend_str = legend_str + newline;
        end
        legend_str = legend_str ...
        + sprintf("$\\quad(R = %.2f, %s$,", correlation_r, p_str) ...
        + newline ...
        + sprintf("$\\quad \\mathrm{NRMSE} = %.1f \\%%, S = %.1f \\%%$)", ...
        100*nrmse, 100*fit_success_ratio);
    end
end

% Plot data points
plot(r_eff_reference', r_eff_median_values, options.marker, ...
    "Color",  options.color, "MarkerEdgeColor", 'k', ...
    "MarkerFaceColor",  options.color, "MarkerSize", options.marker_size, ...
    'DisplayName', legend_str, "HandleVisibility", "on");

% Plot regression line
if options.plot_regression
    regression_str = sprintf("%sLin. Reg.: $%.2f + %.2f \\cdot r_{\\mathrm{eff}}$", ...
        options.legend_prefix, linear_model.Coefficients.Estimate(1), ...
        linear_model.Coefficients.Estimate(2));
    plot(xlims', predict(linear_model, sort(xlims')), ...
        options.regression_line_style, 'Color', options.color, ...
        'LineWidth', 2, 'DisplayName', regression_str);
end

% Plot unity line
if options.plot_unity_line
    existing_unity_line = findobj(gca, 'DisplayName', 'Unity');
    if ~isempty(existing_unity_line)
        delete(existing_unity_line);
    end
    unity_line_lims = [max([xlims(1), ylims(1)]), min([xlims(2), ylims(2)])];
    plot(unity_line_lims, unity_line_lims, 'k--', 'DisplayName', 'Unity', 'HandleVisibility', "on");
end

% Final plot adjustments
pbaspect([xlims(2)-xlims(1) ylims(2)-ylims(1) 1]);
xlabel('$r_{\mathrm{eff}}$ [$\mu$m] (histology)', 'Interpreter', 'latex');
ylabel("$r_{\mathrm{eff}}$ [$\mu$m]", 'Interpreter', 'latex');  
xticks(xtick_vals);
yticks(ytick_vals);
xlim(xlims);
ylim(ylims);
box on;

end
