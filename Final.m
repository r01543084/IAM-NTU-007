%%
%stability 所撰寫之linear advection equation (velocity >0)
%Atmosphere @NTU 2014.1.6
%%
%輸入條件
clear all
close all
clc
x0   = -2*pi        ;%X初始位置
xEnd = 2*pi         ;%X結束位置
dx   = 0.1          ;%每dx切一格
tEnd = 0.1           ;%從0開始計算tEnd秒 
cfl  = 0.04          ;%CFL number
a    = 10            ;%constant
dt   = cfl*dx/a     ;%每dt切一格
scheme_type = 6    ;%算則選擇
inital_type = 1    ;%初始條件，1為step func，2為sin wave
%%
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
%%
switch inital_type
    case 1  %setp func
        u  = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%初始條件
        uu = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%比較用
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
        uu = exp(-(x+4-a*dt*n).^2)+heaviside(x-a*dt*n)-heaviside(x-3-a*dt*n);%比較用
        u(1) = 0;%von Newmann BC
        u(end) = 0;
    case 2  %sin wave
        uu   = sin(x-a*dt*n);%比較用
        u(end-5:end) = uu(end-5:end);%periodic BC
        u(1:5) = uu(1:5);
    end
    
    plot(x,u,'-O',x,uu,'--')
    axis([x0, xEnd, min(uu)-0.5, max(uu)+0.5])
    grid on
    xlabel('x');%水平座標名稱
    ylabel('u(t)');%垂直座標名稱
    title(['time(t) = ',num2str(t(n))]); % 圖形標題
	pause (dt)
end