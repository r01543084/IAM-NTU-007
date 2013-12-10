%%
%Modified Box Methods 1D heat equation
%use Thomas Algorithm solver
%by Atmosphere @ NTU 2013.12.9

%%
clear all;close all;clc;
%時間空間、離散
dx = 1/5;
dt = 0.01;
x0 = 0;
xEnd = 1;
tEnd = 3;
k = 0.1;%heat coef.
x = x0:dx:xEnd;
t = 0:dt:tEnd;
r = k*dt/(2*dx^2);
u = ones(1,length(x));

%Thomas Algorithm solver
dj = ones(1,length(u)-2)*(1+2*r);
bj = ones(1,length(u)-3)*(-r);
aj = ones(1,length(u)-3)*(-r);

%bc
u(1) = -1;%左端點溫度
u(end) = 100;%右端點溫度

%%
%mean loop
for i = 1:length(t)

    
%     A = diag(ones(1,length(x))*(1+2*r),0)+...
%         diag(ones(1,length(x)-1)*(-r),1) +...
%         diag(ones(1,length(x)-1)*(-r),-1);
    
    um  = circshift(u,[0 -1]);%u2
    umm = circshift(u,[0 -2]);%u3
    
    c = r*u+(1-2*r)*um+r*umm;
    c(end-1:end) = [];%最後兩個不可算在內
    c(1) = c(1)+r*u(1);
    c(end) = c(end)+r*u(end);
    
    %===============solving loop==========
    dk = dj(2:end)-bj.*aj./dj(1:end-1);
    dk = [dj(1) dk];
    
    ck = c(2:end) - ( bj.*c(1:end-1)./dk(1:end-1) );
    ck = [c(1) ck];
    
    u(end-2:-1:2) = (ck(end-1:-1:1)-(aj(end:-1:1).*u(end-1:-1:3)))...
                     ./dk(end-1:-1:1);%main equation
                 
	u(end-1) = ck(end)/dk(end);%back substitution according
    
	plot(x,u,'.')
    grid on
    pause(dt)
    
end