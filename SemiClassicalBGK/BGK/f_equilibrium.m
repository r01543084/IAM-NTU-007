function f = f_equilibrium(z,marco_u,mirco_v,T,theta) 
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

f = 1./(exp( (mirco_v-marco_u).^2 ./ T)./z + theta);