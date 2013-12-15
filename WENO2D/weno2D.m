%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%             2D-WENO3   5th order by RK4                                 %
%                      Atmosphere@ntu 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%periodic bc in LF_flux
%輸入條件
clear all;clc;
x0   = -2*pi        ;%X初始位置
xEnd = 2*pi         ;%X結束位置
y0   = -2*pi        ;%y初始位置
yEnd = 2*pi         ;%y結束位置
dx   = 0.3          ;%每dx切一格
dy   = 0.3          ;%每dy切一格
tEnd = 30           ;%從0開始計算tEnd秒
u = 10              ;%x方向速度
v = 10              ;%y方向速度
cflx = 1           ;%x方向CFL number
cfly = 1           ;%y方向CFL number
dt = min(cflx*dx/u,cfly*dy/v);%因u=\=v，會使得cfl 過於太大。因為dt 只有1個，
                              %但是會影響x方向及y方向，故需重新計算dt and cfl
cflx = dt*u/dx      %renew cfl x direction
cfly = dt*v/dy      %renew cfl y direction

%%
x = x0 : dx : xEnd;%切x方向空間網格
y = y0 : dy : yEnd;%切y方向空間網格
t = 0  : dt : tEnd;%切時間網格

%%
%flux type
type = 2;%(1)Linear Advection, (2)Burgers' equation

%initial condition
[X,Y] = meshgrid(x,y);
U = exp(-(X).^2-(Y).^2);
%%
%mean progream
for n = 1:length(t)
    
    for i = 1:length(y)
    %rk method 4th order
    k1 = ( (-LF_flux(type,u,weno3(U(i,:))))        /dx );
    k2 = ( (-LF_flux(type,u,weno3(U(i,:)+k1/2*dt)))/dx );
    k3 = ( (-LF_flux(type,u,weno3(U(i,:)+k2/2*dt)))/dx );
    k4 = ( (-LF_flux(type,u,weno3(U(i,:)+k3*dt)))  /dx );
	U(i,:) = U(i,:) + 1/6*(k1+2*k2+2*k3+k4)*dt;%mean equation
    U(1,:) = (U(end,:)+U(1,:))/2;%periodic bc
    end
    
	for j = 1:length(x)
    %rk method 4th order
    k1 = ( (-LF_flux(type,v,weno3(U(:,j)')))        /dx );
    k2 = ( (-LF_flux(type,v,weno3(U(:,j)'+k1/2*dt)))/dx );
    k3 = ( (-LF_flux(type,v,weno3(U(:,j)'+k2/2*dt)))/dx );
    k4 = ( (-LF_flux(type,v,weno3(U(:,j)'+k3*dt)))  /dx );
	U(:,j) = U(:,j) + 1/6*(k1+2*k2+2*k3+k4)'*dt;%mean equation
    U(:,1) = (U(:,end)+U(:,1))/2;%periodic bc
    end
    contourf(x,y,U)
    grid on
    pause(dt)
end