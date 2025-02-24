function add_subfigure_labels(plot_axes, mode, options)
% add_subfigure_labels Adds subfigure labels (a, b, c...) to a grid of axes.
%
%   This function labels subfigures within a given grid of axes using either 
%   character labels ('a', 'b', 'c'...) or a combination of character and number 
%   labels ('a1', 'a2', 'b1'...).
%
% USAGE:
%   add_subfigure_labels(plot_axes, mode, "position", [x, y], "background", [R, G, B]);
%
% INPUTS:
%   plot_axes - (MxN array of axes handles) Target axes for labeling.
%   mode      - (string) Labeling mode:
%               * "char"    - Uses letters (a, b, c, ...).
%               * "charnum" - Uses letter-number pairs (a1, a2, ...).
%   options:
%     * position    (1x2 vector, default: [0.02, 0.955]) - Label position in normalized coordinates.
%     * background  (1x3 RGB vector, default: []) - Background color for labels.

    arguments
        plot_axes (:,:);
        mode (1,1) string = [];
        options.position = [0.02, 0.955];
        options.background = [];
    end

    charlbl = compose("%s",('a':'z').'); 
    numlbl = compose("%s",(1:9).');

    plot_settings = struct();
    if ~isempty(options.background)
        plot_settings.BackgroundColor = options.background;
    end
    cellargs = namedargs2cell(plot_settings);

    count = 1;
    for y = 1:size(plot_axes,1)
        for x = 1:size(plot_axes,2)     
            
            if mode == "char"
                label_text = charlbl(count);
            elseif mode == "charnum"
                label_text = sprintf("%s%d", charlbl(y), x);
            end
            
            if (plot_axes(y,x).Visible == "off")
                continue;
            end
            count = count + 1;
            text(plot_axes(y,x), ...
                options.position(1), ...
                options.position(2), ...
                label_text, ...
                'Units', 'normalized', ...
                'FontSize', 16, ...
                "FontWeight", "bold", cellargs{:});
        end
    end

end