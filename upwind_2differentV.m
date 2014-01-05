clear all
close all
clc
tEnd = 3         ;%�q0�}�l�p��tEnd�� 
dt   = 0.01     ;%�Cdt��p��@��
x0   = -pi        ;%X��l��m
xEnd = pi         ;%X������m
dx   = 0.01     ;%�Cdx���@��
%C    = abs(u)*dt/dx   ;%CFL number
type = 2         ;%��l����A1��step func�A2��sin wave
x = x0 : dx : xEnd;%���Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����
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