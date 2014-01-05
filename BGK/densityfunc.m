function [n,marco_u1,T] = densityfunc(f,weight,marco_u,mirco_v,Tin,density,dv,nv)
%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)
     % density
     n = sum(weight .* f);%u
     n = repmat(n,nv,1);
     n = sum(weight .* n.*exp(-(0-mirco_v).^2./(2*Tin*0.287)));%v
     n = repmat(n,nv,1);
     n = sum(weight .* n.*exp(-(0-mirco_v).^2./(2*Tin*0.287)));%w
     n = repmat(n,nv,1);
     
     % Macrospic velocity
     j_x = sum(mirco_v .* weight .* f);
     j_x = repmat(j_x,nv,1);
     j_x = sum(weight .* j_x.*exp(-(0-mirco_v).^2./(2*Tin*0.287)));
     j_x = repmat(j_x,nv,1);
     j_x = sum(weight .* j_x.*exp(-(0-mirco_v).^2./(2*Tin*0.287)));
     j_x = repmat(j_x,nv,1);
     marco_u1 = j_x./n;%marco u output
     
     % Temperture
     T1 = sum(weight.*(marco_u-mirco_v).^2.*f);
     T1 = repmat(T1,nv,1);
     T1 = sum(weight.*T1.*exp(-(0-mirco_v).^2./(2*Tin*0.287)));
     T1 = repmat(T1,nv,1);
     T1 = sum(weight.*T1.*exp(-(0-mirco_v).^2./(2*Tin*0.287)));
     T1 = repmat(T1,nv,1);
     
     T2 = sum(weight.*f);
     T2 = repmat(T2,nv,1);
     T2 = sum(weight.*(0-mirco_v).^2.*exp(-(0-mirco_v).^2./(2*Tin*0.287)).*T2);
     T2 = repmat(T2,nv,1);
     T2 = sum(weight.*exp(-(0-mirco_v).^2./(2*Tin*0.287)).*T2);
     T2 = repmat(T2,nv,1);
     
     T3 = sum(weight.*f);
     T3 = repmat(T3,nv,1);
     T3 = sum(weight.*exp(-(0-mirco_v).^2./(2*Tin*0.287)).*T3);
     T3 = repmat(T3,nv,1);
     T3 = sum(weight.*(0-mirco_v).^2.*exp(-(0-mirco_v).^2./(2*Tin*0.287)).*T3);
     T3 = repmat(T3,nv,1);
     T = (T1+T2+T3)./(dv*density*0.287);
     
     
     