function color_order = get_default_color_order()
% get_default_color_order Returns a predefined set of colors as an Nx3 RGB matrix.
%
%   This function provides a default color order for plotting, where each row 
%   represents an RGB triplet with values normalized to the range [0,1].
%
% USAGE:
%   color_order = get_default_color_order();
%
% OUTPUT:
%   color_order - (Nx3 matrix) RGB values for predefined colors, normalized to [0,1].
%
%   The default colors are:
%     * [1.000, 0.400, 0.800] - Pink
%     * [0.000, 0.600, 0.600] - Teal
%     * [1.000, 0.647, 0.251] - Orange
%     * [0.000, 0.400, 0.800] - Blue
%     * [0.800, 0.200, 0.200] - Red

color_order = 1./255.*[...
    255, 102, 204; 
    0, 153, 153; 
    255, 165, 64; 
    0, 102, 204; 
    204, 51, 51];

end
