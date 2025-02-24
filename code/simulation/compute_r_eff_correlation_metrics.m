function [r_eff_median_values,...
    r_eff_conf_values,...
    R_observed,...
    p,...
    linear_model,...
    nrmse,...
    fit_success_ratio] = compute_r_eff_correlation_metrics(r_eff_ref, r_eff_est, options)
%   Computes metrics and correlation between estimated and reference
%   ground truth effective radii.
%
%   This function evaluates various metrics and the correlation between 
%   reference effective radius (r_eff_ref) and estimated effective radius 
%   (r_eff_est). It computes median and confidence intervals, performs 
%   a linear regression, calculates the normalized root-mean-square error 
%   (NRMSE), and computes Pearson's correlation with associated 
%   significance determined using a Monte Carlo permutation test.
%
% USAGE:
%   [r_eff_median_values, r_eff_conf_values, correlation_r, correlation_p, ...
%    linear_model, nrmse, fit_success_ratio] = compute_r_eff_correlation_metrics( ...
%    r_eff_ref, r_eff_est, "rand_stream", rand_stream);
%
% INPUTS:
%   r_eff_ref    - (Nx1 or NxM array) Reference effective radius values.
%   r_eff_est    - (Nx1 array) Estimated effective radius values from L independent
%                  simulations or observations. If L=1, r_eff_est can be Nx1.
%   options:
%     * rand_stream (RandStream, default: global stream) - Custom random number 
%       stream for permutation tests.
%
% OUTPUTS:
%   r_eff_median_values  - (Nx1 vector) Median r_eff estimates across simulations.
%   r_eff_conf_values    - (Nx2 matrix) 2.5% and 97.5% percentiles for r_eff estimates.
%   R_obs                - (scalar) Observed Pearson correlation coefficient.
%   p                    - (scalar) p-value from Monte Carlo permutation test.
%   linear_model         - (LinearModel object) Linear regression model (fitlm).
%   nrmse                - (scalar) Normalized root-mean-square error.
%   fit_success_ratio    - (scalar) Fraction of valid (non-NaN) values in r_eff_est.

    arguments
       r_eff_ref
       r_eff_est
       options.rand_stream = []
    end
    
    if isempty(options.rand_stream)
        rand_stream = RandStream.getGlobalStream();
    else
        rand_stream = options.rand_stream;
    end
    
    % repeat reference data if M > 1
    M = size(r_eff_est, 2);  % Number of independent estimates (simulated MRI runs)   
    if M > 1
        r_eff_ref = repmat(r_eff_ref, 1, M);
    end
    
    % median and confidence interval
    r_eff_median_values = squeeze(median(r_eff_est, 2, 'omitnan'));
    r_eff_conf_values = squeeze(prctile(r_eff_est, [2.5, 97.5], 2));

    invalid_idx = isnan(r_eff_est);
    valid_idx = ~invalid_idx;
    
    % linear regression
    linear_model = fitlm(r_eff_ref(valid_idx), r_eff_est(valid_idx));
    
    % NRMSE
    NRMSECalc = @(pred, ref) sqrt(mean((pred-ref).^2))./mean(ref);
    r_eff_est_imputed = zeros(size(r_eff_est));
    r_eff_est_imputed(valid_idx) = r_eff_est(valid_idx);
    nrmse = NRMSECalc(r_eff_est_imputed(:), r_eff_ref(:));
    
    % fraction of successfull fits
    fit_success_ratio = sum(valid_idx(:))./numel(r_eff_est);
 
    % set number of Monte Carlo permutation iterations
    if M > 1
        K = 1000;   % Reduce K for simulations
    else
        K = 1000000; % Use full K for single experiment
    end
    
    % Observed correlation from pooled data
    R_observed = corr(r_eff_ref(valid_idx), ...
        r_eff_est(valid_idx));  

    % Monte Carlo Permutation to determine p-value
    R_null = zeros(M, K); 
    for m = 1:M
        % Get current simulated MRI estimate
        Y_roi = r_eff_est(:, m);
        valid_idx = find(~isnan(Y_roi));
        
        if isempty(valid_idx)
            R_null(m,:) = 1;
            continue;
        end
        
        r_eff_ref_valid = r_eff_ref(valid_idx);
        Y_fixed = Y_roi(valid_idx);
    
        % Generate K shuffled versions of r_eff_ref_valid
        r_eff_ref_shuffled = zeros(length(r_eff_ref_valid), K);
        
        for k = 1:K
            r_eff_ref_shuffled(:, k) = r_eff_ref_valid(randperm(...
                rand_stream, length(r_eff_ref_valid)));
        end
        
        % Compute correlation for K shuffled sets
        R_null(m, :) = corr(r_eff_ref_shuffled, Y_fixed);

    end
    
    % Compute p-value as the fraction of all M*K permutations 
    % where R_perm >= R_observed
    p = mean(abs(R_null(:)) >= abs(R_observed));

end

