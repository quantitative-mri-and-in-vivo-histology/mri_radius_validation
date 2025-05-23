% Source: https://github.com/NYU-DiffusionMRI/SMI
%
% Modifications: function name changed from "get_even_SH" to 
%   "compute_even_order_spherical_harmonics"
% 
% Licensed under custom license (see licenses/NYU_non_commercial.txt).
% Authors: Santiago Coelho (santiago.coelho@nyulangone.org), 
%           Jelle Veraart, Els Fieremans, Dmitry Novikov
%  Copyright (c) 2022 New York University
%              
% REFERENCES:
% - Coelho, S. et al., 2022. Reproducibility of the Standard Model of diffusion 
%   in white matter on clinical MRI systems. Neuroimage. doi: 10.1016/j.neuroimage.2022.119290.
% - Novikov, D.S. et al., 2018. Rotationally-invariant mapping of scalar and 
%   orientational metrics of neuronal microstructure with diffusion MRI. NeuroImage 174, 518–538.
% - Reisert, M. et al., 2017. Disentangling micro from mesostructure by diffusion 
%   MRI: A Bayesian approach. NeuroImage 147, 964–975.

function Ylm_n = compute_even_order_spherical_harmonics(dirs,Lmax,CS_phase)
    % Ylm_n = get_even_SH(dirs,Lmax,CS_phase)
    %
    % if CS_phase=1, then the definition uses the Condon-Shortley phase factor
    % of (-1)^m. Default is CS_phase=0 (so this factor is ommited)
    %
    % By: Santiago Coelho
    if size(dirs,2)~=3
        dirs=dirs';
    end
    Nmeas=size(dirs,1);
    [PHI,THETA,~]=cart2sph(dirs(:,1),dirs(:,2),dirs(:,3)); THETA=pi/2-THETA;
    l=0:2:Lmax;
    l_all=[];
    m_all=[];
    for ii=1:length(l)
        l_all=[l_all, l(ii)*ones(1,2*l(ii)+1)];
        m_all=[m_all -l(ii):l(ii)];
    end
    K_lm=sqrt((2*l_all+1)./(4*pi) .* factorial(l_all-abs(m_all))./factorial(l_all+abs(m_all)));
    if nargin==2 || isempty(CS_phase) || ~exist('CS_phase','var') || ~CS_phase
        extra_factor=ones(size(K_lm));
        extra_factor(m_all~=0)=sqrt(2);
    else
        extra_factor=ones(size(K_lm));
        extra_factor(m_all~=0)=sqrt(2);
        extra_factor=extra_factor.*(-1).^(m_all);
    end
    P_l_in_cos_theta=zeros(length(l_all),Nmeas);
    phi_term=zeros(length(l_all),Nmeas);
    id_which_pl=zeros(1,length(l_all));
    for ii=1:length(l_all)
        all_Pls=legendre(l_all(ii),cos(THETA));
        P_l_in_cos_theta(ii,:)=all_Pls(abs(m_all(ii))+1,:);
        id_which_pl(ii)=abs(m_all(ii))+1;
        if m_all(ii)>0
            phi_term(ii,:)=cos(m_all(ii)*PHI);
        elseif m_all(ii)==0
            phi_term(ii,:)=1;
        elseif m_all(ii)<0
            phi_term(ii,:)=sin(-m_all(ii)*PHI);
        end
    end
    Y_lm=repmat(extra_factor',1,Nmeas).*repmat(K_lm',1,Nmeas).*phi_term.*P_l_in_cos_theta;
    Ylm_n=Y_lm';
end