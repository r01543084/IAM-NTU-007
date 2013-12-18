function [n,j_x,epsilon] = densityfunc(weight,f0,micro_v)
%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)

     n = sum(weight .* f0);% Number density
     j_x = sum(micro_v .* weight .* f0);% Macrospic moment in x
     epsilon = sum(1/2*( micro_v.^2 ).* weight .* f0);% Energy Density