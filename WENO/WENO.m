%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               WENO3   5th order                                         %
%                      Atmosphere@NTU 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%輸入條件
clear all;close all;clc;
x0   = -2*pi        ;%X初始位置
xEnd = 2*pi         ;%X結束位置
dx   = 0.1          ;%每dx切一格
tEnd = 0.1           ;%從0開始計算tEnd秒 
CFL  = 0.04          ;%CFL number
a    = 10            ;%constant
dt   = CFL*dx/a     ;%每dt切一格
%%
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
%%
%mean progream
%u = heaviside(x)-heaviside(x-3);%initial condition
u = sin(x);
for j = 1:length(t)
    F=[];%淨空F
    dF=[];%same
    f=[];%same
    u_weno=[];%same
    u_next=[];%淨空u_next
    apha  = max(abs(a));%微分flux
    u_expand = [u(end-1:end) u u(1:2)];
    for i = 1:length(u_expand)-4
        u_weno(:,i) = weno3(u_expand(i:i+4));%use weno3 method
    end
    u_weno = [u_weno(1,end) u_weno(1,:);
              u_weno(2,:) u_weno(2,end)];%periodic bc
          
    f = a*u_weno;%flux term
	F = 0.5*((f(2,:)+f(1,:))-apha*(u_weno(2,:)-u_weno(1,:)));%LF flux
	u_next = u - ( F(2:end)-F(1:end-1) )/dx*dt;%mean equation
    %uu = heaviside(x-a*dt*j)-heaviside(x-3-a*dt*j);
    uu = sin(x-a*dt*j);
    u_next(1)=u(end);%bc
    u = u_next;%renew u
    plot(x,u,'.',x,uu,'O')
    axis([x0-1, xEnd+1, min(uu)-0.5, max(uu)+0.5])
    xlabel('x');%水平座標名稱
    ylabel('u(t)');%垂直座標名稱
    title(['time(t) = ',num2str(t(j))]); % 圖形標題
    grid on
    pause(dt)
    %error(j) = sum(abs(uu-u));
end