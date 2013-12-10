%%
%Crank-Nicolson 1D heat equation
%use Thomas Algorithm solver
%by Atmosphere @ NTU 2013.12.9

%%
clear all;close all;clc;
%時間空間、離散
dx = 0.5;
dt = 0.01;
x0 = 0;
xEnd = 10;
tEnd = 10;
k = 20;%heat coef.
x = x0:dx:xEnd;
t = 0:dt:tEnd;
r = k*dt/(2*dx^2)
u = ones(1,length(x));

%Thomas Algorithm solver
dj = ones(1,length(u)-2)*(1+2*r);
bj = ones(1,length(u)-3)*(-r);
aj = ones(1,length(u)-3)*(-r);

dk = dj(2:end)-(bj.*aj./dj(1:end-1));
dk = [dj(1) dk];
    


%bc
u(1) = 1;%左端點溫度
u(end) = 10;%右端點溫度

%%
%mean loop
for i = 1:length(t)
    
%     A = diag(ones(1,length(x)-2)*(1+2*r),0)+...
%         diag(ones(1,length(x)-3)*(-r),1) +...
%         diag(ones(1,length(x)-3)*(-r),-1);

    um  = circshift(u,[0 -1]);%u2
    umm = circshift(u,[0 -2]);%u3
    uT = u(2:end-1);%modified CN
    c = r*u+(1-2*r)*um+r*umm;
    c(end-1:end) = [];%最後兩個不可算在內
    c(1) = c(1)+r*u(1);
    c(end) = c(end)+r*u(end);

    %===============solving loop==========
    
    ck = c(2:end) - ( bj.*c(1:end-1)./dk(1:end-1) );
    ck = [c(1) ck];
    
    uT(end) = ck(end)/dk(end);%back substitution according
    
    uT(end-1:-1:1) = (ck(end-1:-1:1)-(aj(end:-1:1).*uT(end:-1:2)))...
                     ./dk(end-1:-1:1);%main equation
    u = [u(1) uT u(end)];
    
        
%     dk = dj(2:end)-bj.*aj./dj(1:end-1);
%     dk = [dj(1) dk];
%     
%     ck = c(2:end) - ( bj.*c(1:end-1)./dk(1:end-1) );
%     ck = [c(1) ck];
%     
%     u(2:end-2) = (ck(1:end-1)-(aj(1:end).*uT))...
%                      ./dk(1:end-1);%main equation
%                  
% 	u(end-1) = ck(end)/dk(end);%back substitution according

%    u(2:end-1) = thomas(aj,bj,c,dj,uT);
    
	plot(x,u,'-*')
    grid on
    pause(dt)
    
end
    
    