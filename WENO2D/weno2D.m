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
y0 = -2*pi          ;%y初始位置
yEnd = 2*pi         ;%y結束位置
dx   = 0.1          ;%每dx切一格
dy = 0.3            ;%每dy切一格
tEnd = 30           ;%從0開始計算tEnd秒
u = 10              ;%x方向速度
v = 10              ;%y方向速度
cflx = .1           ;%x方向CFL number
cfly = .1           ;%y方向CFL number
dt = min(cflx*dx/u,cfly*dy/v);%因u=\=v，會使得cfl 過於太大。因為dt 只有1個，
                              %但是會影響x方向及y方向，故需重新計算dt and cfl
cflx = dt*u/dx      ;%renew cfl x direction
cfly = dt*v/dy      ;%renew cfl y direction

%%
x = x0 : dx : xEnd;%切x方向空間網格
y = y0 : dy : yEnd;%切y方向空間網格
t = 0  : dt : tEnd;%切時間網格

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
    xlabel('x');%水平座標名稱
    ylabel('u(t)');%垂直座標名稱
    title(['time(t) = ',num2str(t(j))]); % 圖形標題
    grid on
    pause(dt)
end