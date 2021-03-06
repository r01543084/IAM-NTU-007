function F = dfdx(f,mirco_v,idvx,nvx,dx,dt,index)
%this flux function is use weno3 (with 5 order) and using RK4th

        %AP. RK step 1
        k1 = ( (-LF_flux(mirco_v,weno3(f,index),idvx,nvx))        /dx );

        %AP. RK step 2
        k2 = ( (-LF_flux(mirco_v,weno3(f+k1/2*dt,index),idvx,nvx))/dx );

        %AP. RK step 3
        k3 = ( (-LF_flux(mirco_v,weno3(f+k2/2*dt,index),idvx,nvx))/dx );

        %AP. RK step 4
        k4 = ( (-LF_flux(mirco_v,weno3(f+k3*dt,index),idvx,nvx))  /dx );
        
        F = 1/6*(k1+2*k2+2*k3+k4);