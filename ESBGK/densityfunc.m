function [n,j_x,epsilon] = densityfunc(weight,f0,micro_v)
%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)

     n = sum(weight .* f0);% Number density
     j_x = sum(micro_v .* weight .* f0);% Macrospic moment in x
     epsilon = sum(1/2*( micro_v.^2 ).* weight .* f0);% Energy Density