function [r_eff, s_powder_average] = simulate_experimental_r_eff_estimation( ...
    s, ...
    mr_protocol, ...
    d_0, ...
    options)
% Simulates acquisition and estimates r_eff.
%
%   Simulates the experimental process by adding optional Gaussian or Rician noise, 
%   computing the powder-averaged signal, and estimating the effective radius (r_eff). 
%   Supports multi-shell and two-shell acquisition protocols, applying either an 
%   analytical power-law ratio method or non-linear fitting for r_eff estimation.
%
% USAGE:
%   [r_eff, s_powder_avg] = simulate_experimental_r_eff_estimation(s, ...
%       mr_protocol, d_0, "noise_type", "gaussian", "noise_level", 0.02, ...
%       "powder_average", "rician_ml", "n_noise_realizations", 500);
%
% INPUTS:
%   s           - (double) Input signal data.
%   mr_protocol - (MrProtocol) Acquisition parameters.
%   d_0         - (double) Reference diffusivity.
%
% OPTIONAL PARAMETERS (options structure):
%   'r_eff_estimator'      - (string, default="") Method for estimating r_eff.
%   'powder_average'       - (string, default="") Powder averaging method.
%   'noise_type'           - (string, default="") Type of noise ("gaussian" or "rician").
%   'noise_level'          - (double, default=NaN) Noise standard deviation per shell.
%   'noise_level_estimate' - (double, default=NaN) Estimated noise level.
%   'f_im_estimate'        - (double, default=0) Correction term for powder averaging.
%   'n_noise_realizations' - (integer, default=1000) Number of noise realizations.
%   'rand_stream'          - (RandStream, default=[]) Custom random stream.
%
% OUTPUTS:
%   r_eff           - (array) Estimated effective axon radius.
%   s_powder_average - (array) Powder-averaged signal.
%
    arguments
       s double
       mr_protocol MrProtocol
       d_0 double
       options.r_eff_estimator (1,1) string = ""
       options.powder_average (1,1) string = ""
       options.noise_type (1,1) string = ""
       options.noise_level double = nan
       options.noise_level_estimate double = nan
       options.f_im_estimate double = 0
       options.n_noise_realizations = 1000
       options.rand_stream = []
       options.spherical_harmonics_order = 6
    end
    
    if isempty(options.rand_stream)
        rand_stream = RandStream.getGlobalStream();
    else
        rand_stream = options.rand_stream;
    end

    s_original_shape = size(s);    
    s = reshape(s, size(s, 1), []);
    
    % add noise if desired
    if ~isempty(options.noise_type) && options.noise_type ~= ""
        % different noise level per shell
        if length(options.noise_level) > 1
            s_noisy = zeros([size(s), options.n_noise_realizations]);
            for shell_index = 1:mr_protocol.num_shells
                unique_bval_mask = ...
                    mr_protocol.bvec_idx_per_shell{shell_index};
                s_roi = s(unique_bval_mask,:);

                noise_gen = @() rand_stream.randn([size(s_roi), ...
                    options.n_noise_realizations]) ...
                    * options.noise_level(shell_index);
                if options.noise_type == "gaussian"
                    s_roi = s_roi + noise_gen();
                elseif options.noise_type == "rician"
                    s_roi = s_roi + complex(noise_gen(), noise_gen());
                    s_roi = abs(s_roi);
                else
                    error("Unknown noise type.");
                end

                s_noisy(unique_bval_mask,:,:) = s_roi;
            end
            s = s_noisy;            
        else
           % same noise level across shells
            noise_gen = @() rand_stream.randn([size(s), ...
                options.n_noise_realizations])...
                * options.noise_level;

            if options.noise_type == "gaussian"
                s = s + noise_gen();
            elseif options.noise_type == "rician"
                s = s + complex(noise_gen(), noise_gen());
                s = abs(s);
            else
                error("Unknown noise type.");
            end
        end
    else
        % no noise
        % s = repmat(s, ...
        %     [ones(1, ndims(s)), options.n_noise_realizations]);
        options.n_noise_realizations = 1;
    end

    s = reshape(s, size(s, 1), []);
    s_powder_average = zeros(mr_protocol.num_shells, size(s,2));
    for shell_index = 1:mr_protocol.num_shells
        unique_bval_mask = mr_protocol.bvec_idx_per_shell{shell_index};
        s_roi = s(unique_bval_mask,:);

        if options.powder_average == "rician_ml"               
            if length(options.noise_level) > 1
                % varying noise level per shell
                s_powder_average(shell_index,:) = ...
                compute_powder_average(...
                    s_roi, ...
                    mr_protocol.bvec(unique_bval_mask,:), ...
                    options.spherical_harmonics_order, ...
                    'sigma', repmat( ...
                    options.noise_level_estimate(shell_index), ...
                    1, size(s_roi,2)));
            else
                % same noise level across shells
                s_powder_average(shell_index,:) = ...
                    compute_powder_average(...
                    s_roi, ...
                    mr_protocol.bvec(unique_bval_mask,:), ...
                    options.spherical_harmonics_order, ...
                    'sigma',repmat( ...
                    options.noise_level_estimate, ...
                    1, size(s_roi,2)));
            end
        elseif options.powder_average == "gaussian_ml"  
            s_powder_average(shell_index,:) = ...
                compute_powder_average(...
                s_roi, ...
                mr_protocol.bvec(unique_bval_mask,:), ...
                options.spherical_harmonics_order);
        else
            error("Unknown powder-average estimator.");
        end
    end
    
    s_powder_average = s_powder_average - options.f_im_estimate;
    s_powder_average(s_powder_average<0) = 0;

    if mr_protocol.num_shells == 2
        % use analytical expression for estimation
        if numel(unique(mr_protocol.gradient_duration_per_shell)) > 1
            error("Found non-unique gradient duration across shells");
        end
        if numel(unique(mr_protocol.gradient_separation_per_shell)) > 1
            error("Found non-unique gradient separation across shells");
        end
        r_eff = estimate_r_eff_power_law_ratio( ...
            s_powder_average, ...
            mr_protocol.bval_per_shell, ...
            mr_protocol.gradient_duration_per_shell(1), ...
            mr_protocol.gradient_separation_per_shell(1), ...
            d_0);
    else
        % use non-linear fitting
        r_eff = estimate_r_eff_multi_shell( ...
            s_powder_average', ...
            mr_protocol.bval_per_shell, ...
            mr_protocol.gradient_duration_per_shell, ...
            mr_protocol.gradient_separation_per_shell, ...
            d_0);
    end

    r_eff = reshape(r_eff, [s_original_shape(2:end), ...
        options.n_noise_realizations]);
    s_powder_average = reshape(s_powder_average, ...
        [mr_protocol.num_shells s_original_shape(2:end), ...
        options.n_noise_realizations]);
end

