%%
%stability 所撰寫之linear advection equation (velocity >0)
%Atmosphere @NTU 2014.1.6
%%
%輸入條件
clear all
close all
clc
tEnd = 5           ;%從0開始計算tEnd秒 
dt   = 0.001         ;%每dt秒計算一次
x0   = -1         ;%X初始位置
xEnd = 1          ;%X結束位置
dx   = 0.01         ;%每dx切一格
a    = 0.1         ;%constan(波速)
cfl    = a*dt/dx   ;%CFL number
scheme_type = 6    ;%算則選擇
inital_type = 1    ;%初始條件，1為step func，2為sin wave
%%
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
%%
switch inital_type
    case 1  %setp func
        u  = heaviside(x+10)-heaviside(x);%初始條件
        uu = heaviside(x+10)-heaviside(x);%比較用
    case 2  %sin wave
        u  = sin(x);%初始條件
        uu = sin(x);%比較用
end

%%
if scheme_type == 6
    aj = ones(1,length(x)-3)*cfl/4;
    bj = -ones(1,length(x)-3)*cfl/4;
    dj = ones(1,length(x)-2);
else
    aj=0;bj=0;dj=0;%讓a b c d 有值
end

%% Main Loop
for n = 2 : length(t)%離散時間(t)坐標
    
    %mean equation
    u = scheme(u,aj,bj,dj,cfl,scheme_type);
	
    %exact part
    switch inital_type
    case 1  %setp func
        uu = heaviside(x+10-a*dt*n)-heaviside(x-a*dt*n);%比較用
        u(1) = 1;%von Newmann BC
        u(end) = 0;
    case 2  %sin wave
        uu   = sin(x-a*dt*n);%比較用
        u(end-5:end) = uu(end-5:end);%periodic BC
        u(1:5) = uu(1:5);
    end
    
    plot(x,u,'.',x,uu,'--')
    axis([x0, xEnd, min(uu)-0.5, max(uu)+0.5])
    grid on
    xlabel('x');%水平座標名稱
    ylabel('u(t)');%垂直座標名稱
    title(['time(t) = ',num2str(t(n))]); % 圖形標題
    pause (dt)
end