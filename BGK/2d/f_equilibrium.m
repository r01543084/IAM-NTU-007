function [g,h] = f_equilibrium(marco_ux,mirco_vx,marco_uy,mirco_vy,T...
                        ,density,id_space,idvx,idvy) 
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
g = density(id_space)./(2*pi*T(id_space))...
    .*exp(-((mirco_vx(idvx)-marco_ux(id_space)).^2+(mirco_vy(idvy)-marco_uy(id_space)).^2)./(2*T(id_space)));
h = T(id_space).*g;