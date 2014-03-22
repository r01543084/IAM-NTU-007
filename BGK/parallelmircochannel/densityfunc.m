function [n,marco_ux,marco_uy,T,E,p] = densityfunc(g,h,weightx,mirco_vx,weighty,mirco_vy...
                                    ,idvx,idvy)
%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)
     % density
     n_vy = sum(weighty(idvy) .* g,4)         ;%先對vy方向積分
     n = sum(weightx(idvx(:,:,:,1)) .* n_vy,3);%再對vx方向積分得到完整的n(x,y)
     
     % Macrospic velocity
        % Marcospic velocity ux
        jx_vy = sum(weighty(idvy) .* g,4)                     ;%先對vy方向積分
        jx = sum(mirco_vx(idvx(:,:,:,1))...
                          .* weightx(idvx(:,:,:,1)).* jx_vy,3);%再對vx方向積分 
        marco_ux = jx./n                                      ;%marco ux
  

        % Marcospic velocity uy
        jy_vy = sum(mirco_vy(idvy) .* weighty(idvy) .* g,4)   ;%先對vy方向積分
        jy = sum(weightx(idvx(:,:,:,1)).* jy_vy,3)            ;%再對vx方向積分
        marco_uy = jy./n                                      ;%marco uy
        
     %Total Energy
        %part 1:  integral (ux-vx)^2*g
        E1_vy = sum(weighty(idvy) .* g,4)                     ;%先對vy方向積分
        E1 = sum(mirco_vx(idvx(:,:,:,1)).^2 ...
                          .* weightx(idvx(:,:,:,1)).* E1_vy,3);%再對vx方向積分
        %part 2:  integral (uy-vy)^2*g
        E2_vy = sum(mirco_vy(idvy).^2 ...
                    .* weighty(idvy) .* g,4)                  ;%先對vy方向積分
        E2 = sum(weightx(idvx(:,:,:,1)) .* E2_vy,3)           ;%再對vx方向積分
        %part 3:  integral g
        E3_vy = sum(weighty(idvy) .* h,4)                     ;%先對vy方向積分
        E3 = sum(weightx(idvx(:,:,:,1)) .* E3_vy,3)           ;%再對vx方向積分
        %Total it.
        E = (E1+E2+E3)./(n*2);
        
     %Internal Energy
        e = E - 0.5*(marco_ux.^2+marco_uy.^2);
     % Temperture
     T = e*2/3;
     
     %pressure
     p = n.*T;
