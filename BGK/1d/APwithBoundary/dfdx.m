function f = dfdx(f,mirco_v,idvx,nvx,dx,dt,method,dectec)
%this flux function is use weno3 (with 5 order) and using RK4th
switch dectec
    case 'b'% boundary domain
            % RK step 1
            k1 = ( (-LF_flux(mirco_v,intpol(f),idvx,nvx))/dx);

            % RK step 2
            k2 = ( (-LF_flux(mirco_v,intpol([f(:,1:4) f(:,5:6)+k1/2*dt]),idvx,nvx))/dx );

            % RK step 3
            k3 = ( (-LF_flux(mirco_v,intpol([f(:,1:4) f(:,5:6)+k2/2*dt]),idvx,nvx))/dx );

            % RK step 4
            k4 = ( (-LF_flux(mirco_v,intpol([f(:,1:4) f(:,5:6)+k3*dt]),idvx,nvx))  /dx );

            F = 1/6*(k1+2*k2+2*k3+k4);
            
            f(:,5:6) = f(:,5:6) + F*dt;
            f(:,1:4) = [];
    case 'n'% normal domain
    switch method
        case 'RK4'
            % RK step 1
            k1 = ( (-LF_flux(mirco_v,weno3(f),idvx,nvx))        /dx );

            % RK step 2
            k2 = ( (-LF_flux(mirco_v,weno3(f+k1/2*dt),idvx,nvx))/dx );

            % RK step 3
            k3 = ( (-LF_flux(mirco_v,weno3(f+k2/2*dt),idvx,nvx))/dx );

            % RK step 4
            k4 = ( (-LF_flux(mirco_v,weno3(f+k3*dt),idvx,nvx))  /dx );

            F = 1/6*(k1+2*k2+2*k3+k4);
            
            f = f + F*dt;
        case 'SSPRK4'
            % step 1
            k1 = f + 0.391752226571890*dt*(-LF_flux(mirco_v,weno3(f),idvx,nvx))/dx;
            
            % step 2
            k2 = 0.444370493651235*f + 0.555629506348765*k1...
                +0.368410593050371*dt*(-LF_flux(mirco_v,weno3(k1),idvx,nvx))/dx;
            
            % step 3
            k3 = 0.620101851488403*f+ 0.379898148511597*k2...
                +0.251891774271694*dt*(-LF_flux(mirco_v,weno3(k2),idvx,nvx))/dx;
            
            %step 4
            k41= (-LF_flux(mirco_v,weno3(k3),idvx,nvx))/dx;%step 4-1
            k4 = 0.178079954393132*f + 0.821920045606868*k3 + ...%step 4-2
                 0.544974750228521*dt*k41;
             
            % step 5
            f = 0.517231671970585*k2+0.096059710526147*k3...
                +0.063692468666290*dt*k41+0.386708617503269*k4...
                +0.226007483236906*dt*(-LF_flux(mirco_v,weno3(k4),idvx,nvx))/dx;
    end
end
            