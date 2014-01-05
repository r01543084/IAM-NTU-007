function [marco_u,p,energy] = macroproperties1d(n,marco_u,T,nv,dv)
%% Recover Macroscopic Properties
% compute back, fugacity, macroscopic velocities, temperature and pressure.
    % Computing first velocites from the momentum: 
   
% to compute fugacity, temperature and pressure, we need to rely on the
% distribution fucntion that we where are using: MB, FD, BE.

    energy = 1/2*n.*marco_u.^2+dv/2*n.*T*8.314;
    p = n.*T*8.314;
% Using Discrete Ordinate Method:
    energy = repmat(energy,nv,1);
    p = repmat(p,nv,1);
    