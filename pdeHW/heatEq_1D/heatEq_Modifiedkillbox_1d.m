%%
%Modified Box Methods 1D heat equation
%use Thomas Algorithm solver
%by Atmosphere @ NTU 2013.12.9
%%%%%%%%%%%%%%%%%%   WARING ���j���Y���ȯ��for loop!!!  %%%%%%%%%%%%%%%%%%%%%
%%
clear all;close all;clc;
%�ɶ��Ŷ��B����
dx = 0.1;
dt = 0.01;
x0 = 0;
xEnd = 10;
tEnd = 10;
k = 10;%heat coef.
x = x0:dx:xEnd;
t = 0:dt:tEnd;
r = k*dt/(2*dx^2)
u = ones(1,length(x));

%Thomas Algorithm solver
aj = ones(1,length(u)-3)*(dx/dt-2*k/dx);
bj = ones(1,length(u)-3)*(dx/dt-2*k/dx);
dj = ones(1,length(u)-2)*(dx/dt+dx/dt+2*k/dx+2*k/dx);
for i = 2:length(dj)
    dj(i) = dj(i)-bj(i-1)/dj(i-1)*aj(i-1);%get new dj
end

%bc
u(1) = 10;%�����I�ū�
u(end) = 1;%�k���I�ū�

%%
%mean loop

for i = 1:length(t)


    up = circshift(u,[0 -1]);%uj-1
    um = circshift(u,[0 +1]);%uj+1
    uT = u(2:end-1);%modified CN
    
    c = 2*k*(um-u)/dx+2*k*(up-u)/dx+(u+um)*dx/dt+(up+u)*dx/dt;
    c(1) = []; c(end) = [];%�h�Y�h��
    c(1) = c(1)-(dx/dt-2*k/dx)*u(1);%c2-b2*u1(n+1)
    c(end) = c(end)-(dx/dt-2*k/dx)*u(end);%cN-1 - aN-1*uN(n+1)
    
    %===============solving loop==========
    
    u(2:end-1) = thomas(aj,bj,c,dj,uT);
    
	plot(x,u,'-O')
    grid on
	xlabel('x');%�����y�ЦW��
    ylabel('T(t,x)');%�����y�ЦW��
    dt*i
    pause(dt)
    
end