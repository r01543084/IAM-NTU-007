function [u] = scheme(u,cfl,type)
    switch type
        case 1 % first order upwind method
            u(2:end) = u(2:end) - cfl*(u(2:end)-u(1:end-1));
            
        case 2 % MacMormack method
            % a variation of the two step Lax-Wendroff scheme
            
                %first step
                u_temp = u(1:end-1) - cfl*(u(2:end)-u(1:end-1));
                u_temp = [u_temp u(end)];
                %second step
                u(2:end) = 0.5*(u(2:end)+u_temp(2:end)...
                                -cfl*(u_temp(2:end)-u_temp(1:end-1)));
        case 3 % Lux-Wendroff scheme
                u(2:end-1) = u(2:end-1) - cfl/2*(u(3:end)-u(1:end-2))...
                           + cfl^2/2*(u(3:end)-2*u(2:end-1)+u(1:end-2));
        case 4 % Warming-Kutler-Lomax scheme
            w = 4*cfl^2-cfl^4;
            % step 1
            u_temp = u(1:end-1)-2/3*cfl*(u(2:end)-u(1:end-1));
            u_temp = [u_temp u(end)];
            % step 2
            u_temp(2:end) = 0.5*(u(2:end)+u_temp(2:end)...
                                -2/3*cfl*(u_temp(2:end)-u_temp(1:end-1)));
            u_temp(1)=u(1);
            % step 3
            u(3:end-2) = u(3:end-2) ...
                         - cfl/24*(-2*u(5:end)+7*u(4:end-1)...
                                 -7*u(2:end-3)+2*u(1:end-4))...
                         - 3/8*cfl*(u_temp(4:end-1)-u_temp(2:end-3))...
                         - w/24*(u(5:end)-4*u(4:end-1)+6*u(3:end-2)...
                                 -4*u(2:end-3)+u(1:end-4));
        case 5 %Beam-Warming scheme
            u(3:end) = u(3:end)-cfl/2*(3*u(3:end)-4*u(2:end-1)+u(1:end-2))...
                               +cfl.^2/2*(u(3:end)-2*u(2:end-1)+u(1:end-2));
                                    
                                    
                
    end
            