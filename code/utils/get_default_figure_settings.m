function figure_settings = get_default_figure_settings()
% Returns a struct with default figure properties.
%
%   This function provides a struct with predefined settings for figure 
%   appearance, including axes, text, and legend properties.
%
% USAGE:
%   figure_settings = get_default_figure_settings();
%
% OUTPUT:
%   figure_settings - (struct) Struct containing default figure settings:
%                     * DefaultAxesLinewidth (double, default: 0.8) - Line width for axes.
%                     * DefaultAxesVisible (string, default: "on") - Visibility of axes.
%                     * DefaultAxesFontSize (int, default: 10) - Font size for axes labels.
%                     * DefaultTextFontSize (int, default: 10) - Font size for text elements.
%                     * DefaultFigureVisible (string, default: "on") - Figure visibility.
%                     * DefaultLegendFontSize (int, default: 6) - Font size for legend.


    figure_settings = struct();
    figure_settings.DefaultAxesLinewidth =  0.8;
    figure_settings.DefaultAxesVisible = 'on';
    figure_settings.DefaultAxesFontSize = 10;
    figure_settings.DefaultTextFontSize =  10;
    figure_settings.DefaultFigureVisible = 'on';
    figure_settings.DefaultLegendFontSize = 6;
end

