function powder_average = compute_powder_average(signals, dirs, order, options)
% compute_powder_average Computes the powder-averaged signal.
%
%   Estimates the powder-averaged signal using either a Rician or Gaussian 
%   maximum-likelihood estimator. Function fits spherical harmonics (SH) 
%   up to the specified order and iteratively refines the coefficients.
%
%   If 'sigma' is provided in options, the Rician ML estimator is used; 
%   otherwise, the Gaussian ML estimator is applied.
%
% USAGE:
%   powder_avg = compute_powder_average(...
%       signals, dirs, sph_order, struct());   % Gaussian ML
%   powder_avg = compute_powder_average(...
%       signals, dirs, sph_order, struct('sigma', sigma)); % Rician ML
%
% INPUTS:
%   signals - (M,N) double    Measured signal intensities for N voxels
%                             and M diffusion directions.
%   dirs    - (M,3) double    Cartesian coordinates of diffusion directions.
%   order   - (integer)       Maximum SH order for fitting.
%   options - struct          (Optional) Structure containing:
%       * sigma - (1,N) double  Estimated noise standard deviation per direction.
%
% OUTPUTS:
%   powder_average - (1,N) double  Estimated powder-averaged signal per voxel.
%

arguments
   signals (:,:) {mustBeNumeric, mustBeReal}  % Signal matrix
   dirs (:,3) {mustBeNumeric, mustBeReal}     % Diffusion directions
   order (1,1) {mustBeNumeric, mustBeInteger, mustBeReal}  % SH order
   options.sigma (:,:) {mustBeNumeric, mustBeReal} = []  % Noise standard deviation (optional)
end

X = compute_even_order_spherical_harmonics(dirs, order, 0);
% Gaussian ML Estimator 
coef_start = X \ signals;

if isempty(options.sigma)
    % stick with Gaussian ML estimate 
    coef = coef_start;
else
    % Rician ML Estimator adapted from Varadarajan et al.
    % (https://doi.org/10.1109/TMI.2015.2427157)
    sigma = options.sigma;
    sigma(sigma < eps) = eps; 
    signals(signals < eps) = eps;
    max_iter = 50;
    tol = 1e-4;
    coef = zeros(size(coef_start));

    parfor voxel_index = 1:size(signals, 2)  
        W = diag(1 ./ sigma(:, voxel_index).^2 .* ones(size(X, 1), 1));
        P = X' * X * (1 ./ sigma(:, voxel_index).^2);
        R = chol(P);
        neg_XW = -X' * W;
        R_T_inv = R' \ eye(size(R));
        coef_next = coef_start(:, voxel_index);
        coef_next(coef_next < 0) = 0;

        for iter_index = 1:max_iter
            coef_curr = coef_next;
            bessel_arg = (X * coef_curr) .* signals(:, voxel_index) ...
                ./ (sigma(:, voxel_index).^2);
            i0 = besseli(0, bessel_arg);
            i1 = besseli(1, bessel_arg);
            i0(i0 < 0) = eps;
            bessel_ratio = (i1 ./ i0);
            bessel_ratio(isinf(i1) & isinf(i0)) = 1;
            coef_next = R \ (R_T_inv ...
                * (-neg_XW * (signals(:, voxel_index) .* bessel_ratio)));

            if sum(abs(coef_next - coef_curr)) < tol
                break;
            end
        end

        coef(:, voxel_index) = coef_next;
    end
end

powder_average = coef(1, :) ./ sqrt(4 * pi);

end
