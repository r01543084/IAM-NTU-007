function [density2,u2,Ms] = RHcondition(density1,p1,p2,u1,gamma)
% Rankine-Hugoniot condition Reference by 
% http://www.potto.org/gasDynamics/node117.php
    c1 = (gamma*p1/density1);
    Ms = u1+c1*(1+(gamma+1)/(2*gamma)*(p2/p1-1))^0.5
    u2 = u1*( (gamma+1)/(gamma-1)*p2/p1 ) / (1+( (gamma+1)/(gamma-1)*p2/p1 ));
    density2 = density1*u1/u2;
    