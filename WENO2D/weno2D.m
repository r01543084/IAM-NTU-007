%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%             2D-WENO3   5th order by RK4                                 %
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
y0 = -2*pi          ;%y��l��m
yEnd = 2*pi         ;%y������m
dx   = 0.1          ;%�Cdx���@��
dy = 0.3            ;%�Cdy���@��
tEnd = 30           ;%�q0�}�l�p��tEnd��
u = 10              ;%x��V�t��
v = 10              ;%y��V�t��
cflx = .1           ;%x��VCFL number
cfly = .1           ;%y��VCFL number
dt = min(cflx*dx/u,cfly*dy/v);%�]u=\=v�A�|�ϱocfl �L��Ӥj�C�]��dt �u��1�ӡA
                              %���O�|�v�Tx��V��y��V�A�G�ݭ��s�p��dt and cfl
cflx = dt*u/dx      ;%renew cfl x direction
cfly = dt*v/dy      ;%renew cfl y direction

%%
x = x0 : dx : xEnd;%��x��V�Ŷ�����
y = y0 : dy : yEnd;%��y��V�Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����

%%
%flux type
type = 1;%(1)Linear Advection, (2)Burgers' equation

%initial condition
[X,Y] = meshgrid(x,y);
U = exp(-(X).^2-(Y).^2);
%%
%mean progream
for j = 1:length(t)

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
    pause(dt)
end