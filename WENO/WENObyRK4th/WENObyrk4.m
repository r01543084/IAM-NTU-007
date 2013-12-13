%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               WENO3   5th order by RK4                                  %
%                      Atmosphere@ntu 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%periodic bc in LF_flux
%��J����
clear all;clc;
x0   = -2*pi        ;%X��l��m
xEnd = 2*pi         ;%X������m
dx   = 0.1          ;%�Cdx���@��
tEnd = 10           ;%�q0�}�l�p��tEnd�� 
CFL  = 0.5          ;%CFL number
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
%u = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%initial condition
u = sin(x);
for j = 1:length(t)
    F=[];%�b��F
    dF=[];%same
    f=[];%same
    u_weno=[];%same
    u_next=[];%�b��u_next 
    %%
    %rk method 4th order
    k1 = ( (-LF_flux(type,a,weno3(u)))        /dx );
    k2 = ( (-LF_flux(type,a,weno3(u+k1/2*dt)))/dx );
    k3 = ( (-LF_flux(type,a,weno3(u+k2/2*dt)))/dx );
    k4 = ( (-LF_flux(type,a,weno3(u+k3*dt)))  /dx );
	u = u + 1/6*(k1+2*k2+2*k3+k4)*dt;%mean equation
    u(1) = (u(end)+u(1))/2;%periodic bc
    
    %uu = exp(-(x+4-a*dt*j).^2)+heaviside(x-a*dt*j)-heaviside(x-3-a*dt*j);
    uu = sin(x-a*dt*j);
    plot(x,u,'.',x,uu,'O')
    axis([x0-1, xEnd+1, min(uu)-0.5, max(uu)+0.5])
    xlabel('x');%�����y�ЦW��
    ylabel('u(t)');%�����y�ЦW��
    title(['time(t) = ',num2str(t(j))]); % �ϧμ��D
    grid on
    error(j)=sum(abs(uu-u)/length(u));
    pause(dt)
end