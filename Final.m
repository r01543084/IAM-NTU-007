%%
%stability 所撰寫之linear advection equation
%空間離散使用upwind method
%時間離散使用euler forward method
%%
%輸入條件
clear all
close all
clc
tEnd = 3         ;%從0開始計算tEnd秒 
dt   = 0.01     ;%每dt秒計算一次
x0   = -pi        ;%X初始位置
xEnd = pi         ;%X結束位置
dx   = 0.1     ;%每dx切一格
a    = 0.5       ;%constan(波速)
cfl    = a*dt/dx   ;%CFL number
type = 1         ;%初始條件，1為step func，2為sin wave
%%
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
%%
switch type
    case 1  %setp func
        u = heaviside(x);%初始條件
        uu = heaviside(x);%比較用
    case 2  %sin wave
        u = sin(x);%初始條件
        uu   = sin(x);%比較用
end
ap = max(a,0);
am = min(a,0);

for n = 2 : length(t)%離散時間(t)坐標

    u = scheme(u,cfl,1);
	
    %exact part
    switch type
    case 1  %setp func
        uu = heaviside(x-a*dt*n);%比較用
        u(1) = 0;%von Newmann BC
    case 2  %sin wave
        
        uu   = sin(x-a*dt*n);%比較用
        u(1) = u(end);%periodic BC
end
    
    plot(x,u,'*',x,uu,'--')
    axis([x0, xEnd, min(uu), max(uu)])
    grid on
    xlabel('x');%水平座標名稱
    ylabel('u(t)');%垂直座標名稱
    title(['time(t) = ',num2str(t(n))]); % 圖形標題
    pause (dt)
end