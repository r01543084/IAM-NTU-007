%%
%Modified Box Methods 1D heat equation
%use Thomas Algorithm solver
%by Atmosphere @ NTU 2013.12.9

%%
clear all;close all;clc;
%時間空間、離散
dx = 0.2;
dt = 0.01;
x0 = 0;
xEnd = 1;
tEnd = 10;
k = 1;%heat coef.
x = x0:dx:xEnd;
t = 0:dt:tEnd;
u = ones(1,length(x));

%Thomas Algorithm solver
dj = ones(1,length(u)-2)*(dx/dt+dx/dt+2*k/dx+2*k/dx);
bj = ones(1,length(u)-3)*(dx/dt-2*k/dx);
aj = ones(1,length(u)-3)*(dx/dt-2*k/dx);

%bc
u(1) = 1;%左端點溫度
u(end) = 10;%右端點溫度

%%
%mean loop
for i = 1:length(t)

    
%     A = diag(ones(1,length(x))*(1+2*r),0)+...
%         diag(ones(1,length(x)-1)*(-r),1) +...
%         diag(ones(1,length(x)-1)*(-r),-1);

    up = circshift(u,[0 -1]);%uj-1
    um = circshift(u,[0 +1]);%uj+1
    uT = u(2:end-2);%modified CN
    
    c = 2*k*(um-u)/dx+2*k*(up-u)/dx+(u+um)*dx/dt+(up+u)*dx/dt;
    c(1) = []; c(end) = [];%去頭去尾
    c(1) = c(1)-(dx/dt-2*k/dx)*u(1);%c2-b2*u1(n+1)
    c(end) = c(end)-(dx/dt-2*k/dx)*u(end);%cN-1 - aN-1*uN(n+1)
    
    %===============solving loop==========
    
    u(2:end-1) = thomas(aj,bj,c,dj,uT);
    
	plot(x,u,'-O')
    grid on
    pause(dt)
    
end