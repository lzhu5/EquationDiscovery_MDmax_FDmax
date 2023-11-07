function [eta, uvel, theta, zele] = StreamFunction_surface_u (H, h, T, plottype, phase_theta, zlocation)
    %% by Ling Zhu, Feb. 2, 2014

    grav = 9.806 ;

    EorS = 'Eulerian' ;
    uEorS = 0.0 ;
    N = 22 ; 

    NT = 1000 ; 
    NZ = 1000 ; 
    
    [A, B, uBar, cel, wavenumber, R] = ...
           generateStreamFile(H,h,T,EorS,uEorS,N) ; 
    k = wavenumber ; 

    %% surface elevation   
    theta_series = linspace(0, 2*pi, NT) ;
    if (strcmp(plottype, 'temporal')==1 || strcmp(plottype, 'spatial')==1)    
        theta = theta_series ;
    elseif (strcmp(plottype, 'vertical')==1)
        theta = 0 + phase_theta*2*pi ; 
    end

    eta = zeros(size(theta)) + A(end);
    for j=1:N
        eta = eta + A(j) * cos(-j*theta);
    end

    %% velocity: uvel and wvel
    if (strcmp(plottype, 'vertical')==1)        
        zele = linspace(-h, eta, NZ) ;
    else
        id = find(theta>0) ;
        zele = zlocation ;
    end

    j=1:N ;
    if (strcmp(plottype, 'vertical')==1)            
        for iz = 1:length(zele)
            uvel(iz) = (cel-uBar) + ...
                   sum(B(j).*cosh(j*k*(zele(iz)+h)) ./ cosh(j*k*h) .* cos(j*theta));
            wvel(iz) = sum(B(j) .* sinh(j*k*(zele(iz)+h)) ./ cosh(j*k*h) .* sin(j*theta));
        end
    else
        for itheta = 1:length(theta)
            uvel(itheta) = (cel-uBar) + ...
                        sum(B(j).*cosh(j*k*(zele+h)) ./ cosh(j*k*h) .* cos(j*theta(itheta)));
            wvel(itheta) = sum(B(j) .* sinh(j*k*(zele+h)) ./ cosh(j*k*h) .* sin(j*theta(itheta)));        
        end
    end

end