%%
%Crank-Nicolson 1D heat equation
%use Thomas Algorithm solver
%by Atmosphere @ NTU 2013.12.9
%%%%%%%%%%%%%%%%%%   WARING 遞迴關係式僅能用for loop!!!  %%%%%%%%%%%%%%%%%%%%%
%%
clear all;close all;clc;
%時間空間、離散
dx = 0.1;
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
aj = ones(1,length(u)-3)*(-r);
bj = ones(1,length(u)-3)*(-r);
dj = ones(1,length(u)-2)*(1+2*r);

for i = 2:length(u)-2
    dj(i) = dj(i)-bj(i-1)/dj(i-1)*aj(i-1);%get new dj
end

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
    cj = r*u+(1-2*r)*um+r*umm;
    cj(end-1:end) = [];%最後兩個不可算在內
    cj(1) = cj(1)+r*u(1);
    cj(end) = cj(end)+r*u(end);

    %===============solving loop==========

    u(2:end-1) = thomas(aj,bj,cj,dj,uT);
    
	plot(x,u,'-*')
    grid on
    pause(dt)
    
end
    
    