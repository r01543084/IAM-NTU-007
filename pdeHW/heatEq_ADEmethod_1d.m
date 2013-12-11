%%
%ADE method 1D heat equation
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

%Thomas Algorithm solver factor
aj = zeros(1,length(u)-2);
bj = ones(1,length(u)-2)*(-r);
dj = ones(1,length(u)-1)*(1+r);
for i = 2:length(dj)
    dj(i) = dj(i)-bj(i-1)/dj(i-1)*aj(i-1);%get new dj
end

%bc
u(1) = 1;%左端點溫度
u(end) = 10;%右端點溫度

%%
%mean loop

for i = 1:length(t)
    
	um = circshift(u,[0 -1]);%uj-1
    %step 1 start
    p = u;
    uT = p(2:end-1);%modified CN
    
    cp = u*(1-r)+r*um;
    cp(1) = []; cp(end) = [];%去頭去尾
    cp(1) = cp(1)+r*u(1);%c2+b2*u1(n+1)
    
    p(2:end-1) = thomas(aj,bj,cp,dj,uT);
    
    %step 2 start
    q = u;
    uT = q(2:end-1);%modified CN
    
    cq = r*u+(1-r)*um;
    cq(end-1:end) = [];%到N-2個
    cq(end) = cq(end)+r*u(end);%c2+b2*u1(n+1)

    q(2:end-1) = thomas(bj,aj,cq,dj,uT);

    u = (p+q)/2;

	plot(x,u,'-O')
    grid on
    pause(dt)
    
end