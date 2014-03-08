function [n,marco_u,T,E,p] = densityfunc(g,h,weight,mirco_v,idv)
%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)
     % density
     n = sum(weight(idv) .* g);%density

     % Macrospic velocity
     j_x = sum(mirco_v(idv) .* weight(idv) .* g);
     marco_u = j_x./n;%marco u
     
     %Energy
     E = 0.5*(sum(weight(idv).*(mirco_v(idv)).^2.*g)...
                + sum(weight(idv).*h));
     e = E./n-marco_u.^2/2;
     
     % Temperture
     T = e*2/3;
     
     %pressure
     p = n.*T;
