function s_perp = simulate_cylinder_signal_matrix_method(r_bin_centers, ...
    fiber_dir, g, bvec, small_delta, big_delta, d_0, tau)
% Simulates dMRI signal using matrix method.
%
%   Simulates the perpendicular diffusion-weighted MRI signal inside a cylinder 
%   for a PGSE sequence with a rectangular gradient waveform, constant 
%   gradient duration and separation, but variable diffusion gradient amplitude.
%
%   This function is a reimplementation of the Cylinder_GEN.m function from the 
%   MISST Toolbox (Drobnjak et al. 2010, Drobnjak et al. 2011, Ianus et al. 2012, 
%   Ianus et al. 2016, http://mig.cs.ucl.ac.uk/index.php?n=Tutorial.MISST). 
%   It is specialized for PGSE sequences with a rectangular gradient waveform 
%   and constant gradient duration and separation, but variable diffusion 
%   gradient amplitude, for multiple gradient directions and cylinder radii. 
%   Parts copied/modified from the original implementation are marked in the code.
%
% USAGE:
%   s_perp = simulate_cylinder_signal_matrix_method(r_bin_centers, fiber_dir, g, ...
%       bvec, small_delta, big_delta, d_0, tau);
%
% INPUTS:
%   r_bin_centers - (:,1) double  Cylinder radii.
%   fiber_dir     - (:,3) double  Fiber direction unit vectors.
%   g             - (:,1) double  Diffusion gradient amplitudes.
%   bvec          - (:,3) double  Diffusion gradient directions.
%   small_delta   - (double)      Gradient pulse duration.
%   big_delta     - (double)      Gradient separation time.
%   d_0           - (double)      Free diffusivity.
%   tau           - (double)      Time step.
%
% OUTPUT:
%   s_perp - (array) Simulated perpendicular diffusion signal.
%


    % Compute angle of fiber with z-direction (theta) and azimuthal angle (phi).
    phis = zeros(size(fiber_dir,1),1);
    thetas = zeros(size(fiber_dir,1),1);
    for fiber_dir_index = 1:size(fiber_dir,1)
        phi = cart2sph(fiber_dir(fiber_dir_index,1), fiber_dir(fiber_dir_index,2), fiber_dir(fiber_dir_index,3));
        theta = compute_angles(fiber_dir(fiber_dir_index,:), [0 0 1]);
        phis(fiber_dir_index) = phi;
        thetas(fiber_dir_index) = theta;
    end

    %% 
    % Initialize constants with default values
    % Default values as in MISST Toolbox (MMConstants.m)
    gstep = 0.04/100;
    GAMMA = 2.675987E8;
    dimA = 20;
    qunit = GAMMA*tau*gstep;
    qunit = qunit/(2*pi);
   
    % Compute A, S, R, and kDN for cylindrical geometry
    % See MISST Toolbox, reimplemented from MMConstants.m
    N_r = length(r_bin_centers);
    A_cyl = zeros(N_r, dimA, dimA);
    R_cyl = zeros(N_r, dimA, dimA);
    kDn_cyl = zeros(N_r, dimA, dimA);
    S_cyl = zeros(N_r, dimA);
    cylinder_modelcode = 1;
    for j = 1:length(r_bin_centers)
        roo = MMroots(dimA,r_bin_centers(j),cylinder_modelcode); % cylinder
        [S_cyl(j,:), A_cyl(j,:,:), kDn_cyl(j,:,:)] = MatrixSA(qunit,r_bin_centers(j),dimA,roo,cylinder_modelcode);
        R_cyl(j,:,:) = MatrixR(tau,r_bin_centers(j),dimA,roo,cylinder_modelcode,d_0);
    end
    
    % Compute v and w for all cylinder directions.
    % See MISST Toolbox, reimplemented from Cylinder_GEN.m
    N_fibres = length(thetas);
    vs = zeros(length(thetas), 3);
    ws = zeros(length(thetas), 3);
    for i = 1:length(thetas)       
        cos_theta = cos(thetas(i));
        sin_theta = sqrt(1-cos_theta^2);
        cos_phi = cos(phis(i));
        sin_phi = sqrt(1-cos_phi^2);        
        vs(i,:) = [cos_theta*cos_phi^2+sin_phi^2; -(1-cos_theta)*sin_phi*cos_phi; -sin_theta*cos_phi]; % vectors for the new x and y directions; see Ozarslan 2010
        ws(i,:) = [ -(1-cos_theta)*sin_phi*cos_phi; cos_theta*sin_phi^2+cos_phi^2; -sin_theta*sin_phi];
    end
    
    % Determine uniuqe v and w
    [unique_rows, unique_row_ids, row_labels] = unique([ws vs], 'rows');
    N_unique_fibres = length(unique_row_ids);
    vs_unique = unique_rows(:,1:3);
    ws_unique = unique_rows(:,4:6);

    N_bvec = size(bvec,1);
    bvec_dot_ws = bvec*vs_unique';
    bvec_dot_vs = bvec*ws_unique';
    s_perp_unique = zeros(N_r, N_unique_fibres, size(bvec,1), "like", 1j);

    % See Callaghan et al. 1997, Eq. 26 and Fig. 3b
    M = round(small_delta/tau);
    N = round(big_delta/tau);

    % Compute signals per fiber direction and radius
    parfor fiber_index = 1:N_unique_fibres

        for r_index = 1:N_r
            
            R = squeeze(R_cyl(r_index,:,:));
            A_base = squeeze(A_cyl(r_index,:,:));
            S = squeeze(S_cyl(r_index,:,:));
            kDn = squeeze(kDn_cyl(r_index,:,:));    
            R_pow = R^(N-M);  
            
            for m=1:N_bvec
                
                % See MISST toolbox, reimplemented from Cylinder_GEN.m
                A = (real(A_base)+(imag(A_base).*kDn)*bvec_dot_ws(m,fiber_index))+1i*imag(A_base)*bvec_dot_vs(m,fiber_index);
                A = A^(round(g(m)./gstep)); 

                % See Callaghan et al. 1997, Eq. 26
                s_perp_unique(r_index,fiber_index,m) = S*(R*A)^M*R_pow*(R*A')^M*R*S';
            end       
        end
    end  
        
    % Fill duplicate values from unique ones
    s_perp = zeros(N_r, N_fibres, size(bvec,1), "like", 1j);
    for fiber_index = 1:N_fibres
        s_perp(:,fiber_index,:) = s_perp_unique(:,row_labels(fiber_index),:);
    end
   
end