%THOMAS ALGORITHM
function [uk] = thomas(aj,bj,c,dj,u)

    dk = dj(2:end)-bj.*aj./dj(1:end-1);
    dk = [dj(1) dk];
    
    ck = c(2:end) - ( bj.*c(1:end-1)./dk(1:end-1) );
    ck = [c(1) ck];
    
    uk = (ck(1:end-1)-(aj(1:end).*u))...
                        ./dk(1:end-1);%main equation
                 
	uk = [uk ck(end)/dk(end)];%back substitution according
