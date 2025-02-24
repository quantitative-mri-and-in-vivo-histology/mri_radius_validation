function [ax_inset] = plot_inset(ax_full, inset_rel_pos, inset_xlim, inset_ylim)
% Creates an inset axis with a magnified excerpt within an existing plot.
%
%   This function generates an inset within an existing axis, copying its contents 
%   and applying zoomed-in axis limits. It also adds a rectangle annotation to 
%   indicate the inset area in the main plot.
%
% USAGE:
%   ax_inset = plot_inset(ax_full, inset_rel_pos, inset_xlim, inset_ylim);
%
% INPUTS:
%   ax_full       - (axes handle) Parent axis where the inset will be placed.
%   inset_rel_pos - (1x4 vector) Relative inset position [x, y, width, height] 
%                   within the parent axis, in normalized coordinates.
%   inset_xlim    - (1x2 vector) X-axis limits for the inset.
%   inset_ylim    - (1x2 vector) Y-axis limits for the inset.
%
% OUTPUT:
%   ax_inset - (axes handle) Handle to the created inset axis.
%

    ax_full_pos = tightPosition(ax_full);
    inset_pos = [ax_full_pos(1) + ax_full_pos(3) * inset_rel_pos(1), ...
                 ax_full_pos(2) + ax_full_pos(4) * inset_rel_pos(2), ...
                 ax_full_pos(3) * inset_rel_pos(3), ...
                 ax_full_pos(4) * inset_rel_pos(4)];
             
    ax_inset = axes('Position', inset_pos);
    axes(ax_inset);
    copyobj(ax_full.Children, ax_inset);
    xlim(inset_xlim);
    ylim(inset_ylim);
    
    annotation("rectangle", inset_pos, 'Color', "k", 'LineWidth', 1);
    
    axes(ax_full);
    rectangle('Position', [inset_xlim(1), inset_ylim(1), ...
                           inset_xlim(2) - inset_xlim(1), inset_ylim(2) - inset_ylim(1)], ...
              'EdgeColor', "k", 'LineWidth', 1);
    
    axes(ax_inset);
end
