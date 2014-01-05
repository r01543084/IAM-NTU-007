clear all
close all
clc
tEnd = 3         ;%從0開始計算tEnd秒 
dt   = 0.01     ;%每dt秒計算一次
x0   = -pi        ;%X初始位置
xEnd = pi         ;%X結束位置
dx   = 0.01     ;%每dx切一格
%C    = abs(u)*dt/dx   ;%CFL number
type = 2         ;%初始條件，1為step func，2為sin wave
x = x0 : dx : xEnd;%切空間網格
t = 0  : dt : tEnd;%切時間網格
a = zeros(1,length(x));
n=ceil(length(x)/2);
a(1:n)=1;
a(n+1:end)=-1;
ap = max(a,0);%u+
am = min(a,0);%u-
u = exp(-16*(x-2).^2)+exp(-16*(x+2).^2);
for i = t
u(2:end-1) = u(2:end-1) - dt/dx*ap(1:end-2).*(u(2:end-1)-u(1:end-2))...
                        - dt/dx*am(2:end-1).*(u(3:end)-u(2:end-1));
                            plot(x,u,'.')
                            pause(dt)
end