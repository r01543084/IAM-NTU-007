function [n,marco_u,T,E,p] = densityfunc(g,h,weight,mirco_v,idv)
%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)
     % density
     n = sum(weight(idv) .* g);%density

     % Macrospic velocity
     j_x = sum(mirco_v(idv) .* weight(idv) .* g);
     marco_u = j_x./n;%marco u
     
     %Energy
     E = 0.5*(sum(weight(idv).*mirco_v(idv).^2.*g)...
                + sum(weight(idv).*h));
     e = E./n-marco_u.^2/2;
     
     % Temperture
     T = e*2/3;
     
     %pressure
     p = n.*T;
