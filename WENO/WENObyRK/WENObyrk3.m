%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               WENO3   5th order by RK3                                  %
%                      Atmosphere@ntu 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%輸入條件
clear all;close all;clc;
x0   = -2*pi        ;%X初始位置
xEnd = 2*pi         ;%X結束位置
dx   = 0.1          ;%每dx切一格
tEnd = 10           ;%從0開始計算tEnd秒 
CFL  = 0.8          ;%CFL number
a    = 1            ;%constant
dt   = CFL*dx/a     ;%每dt切一格
%%
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
%%
%flux type
type = 1;%(1)Linear Advection, (2)Burgers' equation

%initial condition
%u = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%initial condition
u = sin(x);
%%
%mean progream
for j = 1:length(t)

    %rk method 4th order
    k1 = ( (-LF_flux(type,a,weno3(u)))                 /dx );
    k2 = ( (-LF_flux(type,a,weno3(u+dt*k1/2)))         /dx );
    k3 = ( (-LF_flux(type,a,weno3(u+(k1+k2)*(dt/2)/2)))/dx );
	u = u + 1/6*(k1+k2+4*k3)*dt;%mean equation
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