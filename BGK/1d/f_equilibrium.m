function [g,h] = f_equilibrium(marco_u,mirco_v,T,density,idx,idv) 
% Compute equilibrum in 1d cases for:
%
% MB:  Maxwell-Boltzmann, theta =  0
% FD:  Fermi-Diract,      theta = +1
% BE:  Bose-Einstein,     theta = -1
%
% inputs:
% u: macroscopic velocity
% x: microscopic velocity
% t: temperature
% r: fugacity
%
g = density(idx).*exp(-(mirco_v(idv)-marco_u(idx)).^2./(2*T(idx)))./sqrt(2*pi*T(idx));
%g = density(idx)./(pi*T(idx)).^0.5.*exp(-(marco_u(idx)-mirco_v(idv)).^2./T(idx));
h = 2*T(idx).*g;