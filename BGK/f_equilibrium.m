function f = f_equilibrium(marco_u,mirco_v,T,density,dv) 
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

f = (density./((2*pi*T*0.287).^(dv/2))) .* exp(-((marco_u-mirco_v).^2)./(2*T*0.287));