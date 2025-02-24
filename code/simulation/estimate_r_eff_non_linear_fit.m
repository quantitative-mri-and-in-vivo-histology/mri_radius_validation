% Source: https://github.com/NYU-DiffusionMRI/AxonRadiusMapping
% Original function: "getAxonRadius"
%
% Modifications:
% - Function name changed from "getAxonRadius" to "estimate_r_eff_non_linear_fit"
% - Fixed typo in "Neumann" (now "Neuman")
% - added D0 as parameter
%
% Licensed under the Mozilla Public License 2.0 (see licenses/MPL2.txt for details).
function [r, beta]  = estimate_r_eff_non_linear_fit(delta, Delta, g, y, model, D0)


    % "getAxonRadius": Estimate the axon radius from multi-shell diffusion
    % data, of which the lowest b-value is sufficiently high to suppress
    % the extra-axonal compartment. For in vivo human MRI, we recommend
    % b=6000s/mm2. 
    %
    % [r, beta]  = getAxonRadius(delta, Delta, g, y, model)
    %       output:
    %           - r:  effective MR radius
    %           - beta: [1x1] The slope of the signal scaling as a function
    %           of 1/sqrt. ~ f/sqrt(Da)
    
    %       input:
    %           - delta:  gradient duration [ms]
    %           - Delta:  gradient seperation [ms]
    %           - g:      gradient strenght [mT/m]
    %           - y:      powder-averaged diffusion-weighted signal,
    %                     normalized to y(b=0) = 1. If the dot comparment
    %                     cannot be ignored, it need to be subtracted from y prior to fitting.  
    %           - model:  Van Gelderen or Neuman (the latter only applied in long pulse regimes, but typically applies)
    %
    %  Author: Jelle Veraart (jelle.veraart@nyulangone.org)
    %  Copyright (c) 2019 New York University



    if ~exist('model', 'var')
        model = 'VanGelderen';
    end
        
    
    options = optimset('lsqnonlin'); 
    options = optimset(options,'Jacobian','on','TolFun',1e-12,'TolX',1e-12,'MaxIter',10000,'Display','off');
    
    
    gyroMagnRatio =  267.513*10^(-6);
    q = g*gyroMagnRatio;
    b = (q.*delta).^2.*(Delta - delta/3);
    
    
  
    start = [sqrt(4*pi)*rand(), 1.5+2*rand()];
    pars = lsqnonlin(@(x)residuals(x, [delta(:), Delta(:), g(:)], y, model, D0),start,[0 0],[sqrt(4*pi) 5],options);
    
    beta = pars(1);
    r = pars(2);

end

function [E, J] = residuals(pars, X, y, model, D0)
    delta = X(:, 1);
    Delta = X(:, 2);
    g = X(:, 3);
    
    [shat, dshat] = AxonDiameterFWD(delta, Delta, g, pars, model, D0); 
    E = shat(:) - y(:);
    J = dshat;
    
end

function [s, ds] = AxonDiameterFWD(delta, Delta, g, pars, model, D0)
    
    % read parameters with fixed D0
    gyroMagnRatio =  267.513*10^(-6);
    beta = pars(1);
    r = pars(2);

    % compute b-values using Tjeskall-Tanner sequence
    q = g*gyroMagnRatio;
    b = (q.*delta).^2.*(Delta - delta/3);
 

    % Forward model
    switch model
        case 'Neuman'
            s = beta*exp(-(7/48)*q.^2.*delta*r^4/D0)./sqrt(b);
            ds_dbeta = exp(-(7/48)*q.^2.*delta*r^4/D0) ./ sqrt(b);    
            ds_dr = beta ./ sqrt(b) .* exp(-(7/48)*q.^2.*delta*r^4/D0) .*(-(7/12)*q.^2.*delta*r^3/D0);
            
        case 'VanGelderen'
            [Svg, dSvg] = vg(delta, Delta, q, r, D0);
            s = beta*exp(-Svg) ./ sqrt(b);
            ds_dbeta = exp(-Svg) ./ sqrt(b);
            ds_dr = beta*exp(-Svg) ./ sqrt(b) .* -dSvg;
    end
    
    ds = [ds_dbeta, ds_dr];
end


function [s, ds] = vg(delta, Delta, q, r, D0)

        td=r^2/D0;
        bardelta=delta/td; dbardelta = -2*delta*D0 / r^3;
        barDelta=Delta/td; dbarDelta = -2*Delta*D0 / r^3;


        N=15; 
        b = [1.8412    5.3314    8.5363   11.7060   14.8636   18.0155   21.1644 24.3113   27.4571   30.6019 ...
             33.7462   36.8900   40.0334   43.1766   46.3196 49.4624   52.6050   55.7476   58.8900   62.0323];
   
   
        s = 0; ds = 0;
        for k=1:N
           s = s + (2/(b(k)^6*(b(k)^2-1)))*(-2 + 2*b(k)^2*bardelta + ...
                          2*(exp(-b(k)^2*bardelta)+exp(-b(k)^2*barDelta)) - ...
                          exp(-b(k)^2*(bardelta+barDelta)) - exp(-b(k)^2*(barDelta-bardelta))); 
           ds = ds +  (2/(b(k)^6*(b(k)^2-1)))*( ...
                            2*b(k)^2*dbardelta + ...
                            2*exp(-b(k)^2*bardelta) - b(k)^2*dbardelta + ...
                            2*exp(-b(k)^2*barDelta) - b(k)^2*dbarDelta - ...
                            exp(-b(k)^2*(bardelta+barDelta)) - b(k)^2*(dbardelta + dbarDelta) - ... 
                            exp(-b(k)^2*(barDelta-bardelta)) - b(k)^2*(dbarDelta - dbardelta));           
        end
        s = s.*D0.*q.^2.*td^3;
        ds = ds.*D0.*q.^2.*td^3 + 6*s.*q.^2.*r^5 / D0^2;
end