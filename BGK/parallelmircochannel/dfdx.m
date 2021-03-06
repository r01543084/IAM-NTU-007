function F = dfdx(f,mirco_v,idvx,nvx,dx,dt,index,dector)
%this flux function is use weno3 (with 5 order) and using RK4th
% index for switch x dir. and y dir.
switch dector
    case 'b'
	%AP. RK step 1
        k1 = ( (-LF_flux(mirco_v,intpol(f),idvx,nvx))        /dx );

        %AP. RK step 2
        k2 = ( (-LF_flux(mirco_v,intpol([f(1:4,:,:,:);f(5:6,:,:,:)+k1/2*dt]),idvx,nvx))/dx );

        %AP. RK step 3
        k3 = ( (-LF_flux(mirco_v,intpol([f(1:4,:,:,:);f(5:6,:,:,:)+k2/2*dt]),idvx,nvx))/dx );

        %AP. RK step 4
        k4 = ( (-LF_flux(mirco_v,intpol([f(1:4,:,:,:);f(5:6,:,:,:)+k3*dt]),idvx,nvx))  /dx );
        
        F = 1/6*(k1+2*k2+2*k3+k4);
    case 'n'
        %AP. RK step 1
        k1 = ( (-LF_flux(mirco_v,weno3(f,index),idvx,nvx))        /dx );

        %AP. RK step 2
        k2 = ( (-LF_flux(mirco_v,weno3(f+k1/2*dt,index),idvx,nvx))/dx );

        %AP. RK step 3
        k3 = ( (-LF_flux(mirco_v,weno3(f+k2/2*dt,index),idvx,nvx))/dx );

        %AP. RK step 4
        k4 = ( (-LF_flux(mirco_v,weno3(f+k3*dt,index),idvx,nvx))  /dx );
        
        F = 1/6*(k1+2*k2+2*k3+k4);
end