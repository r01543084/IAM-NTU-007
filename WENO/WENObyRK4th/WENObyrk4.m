%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               WENO3   5th order by RK4                                  %
%                      Atmosphere@ntu 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%��J����
clear all;close all;clc;
x0   = -3*pi        ;%X��l��m
xEnd = 3*pi         ;%X������m
dx   = 0.05          ;%�Cdx���@��
tEnd = 10           ;%�q0�}�l�p��tEnd�� 
CFL  = 0.1          ;%CFL number
a    = 1            ;%constant
dt   = CFL*dx/a     ;%�Cdt���@��
%%
x = x0 : dx : xEnd;%���Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����
%%
%flux type
type = 1;%(1)Linear Advection, (2)Burgers' equation
%%
%mean progream
u = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%initial condition
%u = sin(x);
for j = 1:length(t)
    F=[];%�b��F
    dF=[];%same
    f=[];%same
    u_weno=[];%same
    u_next=[];%�b��u_next 
    %%
    %rk method 4th order
    k1 = dt*( (-LF_flux(type,a,weno3(u)))/dx );
    k2 = dt*( (-LF_flux(type,a,weno3(u+k1/2)))/dx );
    k3 = dt*( (-LF_flux(type,a,weno3(u+k2/2)))/dx );
    k4 = dt*( (-LF_flux(type,a,weno3(u+k3)))/dx );
	u_next = u + 1/6*(k1+2*k2+2*k3+k4);%mean equation
    
    uu = exp(-(x+4-a*dt*j).^2)+heaviside(x-a*dt*j)-heaviside(x-3-a*dt*j);
    %uu = sin(x-a*dt*j);
    u_next(1)=u(end);%bc
    u = u_next;%renew u
    plot(x,u,'.',x,uu,'--')
    axis([x0-1, xEnd+1, min(uu)-0.5, max(uu)+0.5])
    xlabel('x');%�����y�ЦW��
    ylabel('u(t)');%�����y�ЦW��
    title(['time(t) = ',num2str(t(j))]); % �ϧμ��D
    grid on
    %error(j)=sum(abs(uu-u));
    pause(dt)
end