function [z,marco_u,T,pressure] = macroproperties1d(n,j_x,epsilon,nx,nv,theta)
%% Recover Macroscopic Properties
% compute back, fugacity, macroscopic velocities, temperature and pressure.
    % Computing first velocites from the momentum:
    marco_u = j_x./n; 
   
% to compute fugacity, temperature and pressure, we need to rely on the
% distribution fucntion that we where are using: MB, FD, BE.

switch theta
    case{-1} % BE
    % If BE: we apply bisection method to the approx BE distribution Eq.
        r_a = 0.001; r_b = 0.99; tol = 1e-7;L   = 1:50;
        for i = 1:nx
        psi = @(r_x) 2*epsilon(i)-sum((r_x.^L) ./ (L.^1.5))* ...
        (n(i)/sum((r_x.^L) ./ (L.^0.5)))^3/(2*pi)- n(i)*(marco_u(i)^2);
        r_p(i) = bisection(psi,r_a,r_b,tol);
        end
        z = r_p;
        T = n.^2./(pi*(BE(r_p,0.5)).^2);
        pressure = epsilon - 1/2*n.*(marco_u.^2);
        
    case{1} % FD
    % if FD: we apply bisection method to the approx FD distribution Eq.
        r_a = 0.001; r_b = 0.99; tol = 1e-7;
        for i = 1:nx
        psi = @(r_x) 2*epsilon(i)- FD(r_x,1.5)*(n(i)/FD(r_x,0.5))^3/(2*pi) ...
        - n(i)*(marco_u(i)^2);
        r_p = bisection(psi,r_a,r_b,tol);
        z(i) = r_p;
        T(i) = n(i)^2/(pi*(FD(r_p,0.5))^2);
        pressure(i) = epsilon(i) - 1/2*n(i)*(marco_u(i)^2);
        end        
    
    case{0} % MB
    % IF MB: the task is much simple.
        T = 4*epsilon./n - 2*marco_u.^2;
        z = n./sqrt(pi.*T);
        pressure = epsilon - 1/2.*n.*marco_u.^2;
    otherwise 
        error('theta can only be: -1, 0, +1 ');
end

% Using Discrete Ordinate Method:
    z = repmat(z,nv,1); marco_u = repmat(marco_u,nv,1); 
    T = repmat(T,nv,1); %p = repmat(p,nv,1);
    