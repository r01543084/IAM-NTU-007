%%
%stability 所撰寫之linear advection equationGodunov’s method
%空間離散使用upwind method
%時間離散使用euler forward method
%%
%輸入條件
clear all
clc
tEnd = 1       ;%從0開始計算10秒 
dt   = 0.0001    ;%每0.01秒計算一次
x0   = -1       ;%X初始位置
xEnd =  5      ;%X結束位置
dx   = 0.001    ;%每0.1切一格
%%
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
u = exp(-16*x.^2)+1;%初始條件
uCompare = exp(-16*x.^2)+1;%比較用
%%
%主程式
for j = 1:length(t)
    u_next = zeros(size(u));%宣告零變量以免回圈出錯
    for i = 2 : length(x)-1
        u_next(i) = u(i) - dt/dx * (F(u(i),u(i+1))-F(u(i-1),u(i)));%mean equation
    end
    u = u_next;%update info
    u(1) = 1;%BC
    plot(x,u,'*',x,uCompare,'--')
    axis([x0, xEnd, min(uCompare), max(uCompare)])
    grid on
    xlabel('x');%水平座標名稱
    ylabel('u(t)');%垂直座標名稱
    title(['time(t) = ',num2str(t(j))]); % 圖形標題
    pause(dt)
end