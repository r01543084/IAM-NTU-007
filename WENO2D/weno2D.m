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
y0   = -2*pi        ;%y��l��m
yEnd = 2*pi         ;%y������m
dx   = 0.1          ;%�Cdx���@��
dy   = 0.1          ;%�Cdy���@��
tEnd = 30           ;%�q0�}�l�p��tEnd��
u = 10              ;%x��V�t��
v = 1               ;%y��V�t��
cflx = 1            ;%x��VCFL number
cfly = 1            ;%y��VCFL number
dt = min(cflx*dx/u,cfly*dy/v);%�]u=\=v�A�|�ϱocfl �L��Ӥj�C�]��dt �u��1�ӡA
                              %���O�|�v�Tx��V��y��V�A�G�ݭ��s�p��dt and cfl
cflx = dt*u/dx      %renew cfl x direction
cfly = dt*v/dy      %renew cfl y direction

%%
x = x0 : dx : xEnd;%��x��V�Ŷ�����
y = y0 : dy : yEnd;%��y��V�Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����

%%
%flux type
type =1;%(1)Linear Advection, (2)Burgers' equation

%initial condition
[X,Y] = meshgrid(x,y);
U = exp(-(X).^2-(Y).^2)+0.5;
%U = heaviside(X);
%%
%mean progream
for n = 1:length(t)
    
    for i = 1:length(y)
    %rk method 4th order
    k1x = ( (-LF_flux(type,u,weno3(U(i,:))))        /dx );
    k2x = ( (-LF_flux(type,u,weno3(U(i,:)+k1x/2*dt)))/dx );
    k3x = ( (-LF_flux(type,u,weno3(U(i,:)+k2x/2*dt)))/dx );
    k4x = ( (-LF_flux(type,u,weno3(U(i,:)+k3x*dt)))  /dx );
	U(i,:) = U(i,:) + 1/6*(k1x+2*k2x+2*k3x+k4x)*dt;%mean equation
    U(1,:) = (U(end,:)+U(1,:))/2;%periodic bc
    end
    
	for j = 1:length(x)
    %rk method 4th order
    k1y = ( (-LF_flux(type,v,weno3(U(:,j)')))        /dx );
    k2y = ( (-LF_flux(type,v,weno3(U(:,j)'+k1y/2*dt)))/dx );
    k3y = ( (-LF_flux(type,v,weno3(U(:,j)'+k2y/2*dt)))/dx );
    k4y = ( (-LF_flux(type,v,weno3(U(:,j)'+k3y*dt)))  /dx );
	U(:,j) = U(:,j) + 1/6*(k1y+2*k2y+2*k3y+k4y)'*dt;%mean equation
    U(:,1) = (U(:,end)+U(:,1))/2;%periodic bc
    end
    contourf(x,y,U)
    %view(135,45)
    axis tight
    grid on
    pause(dt)
end