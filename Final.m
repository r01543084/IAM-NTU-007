%%
%stability �Ҽ��g��linear advection equation
%�Ŷ������ϥ�upwind method
%�ɶ������ϥ�euler forward method
%%
%��J����
clear all
close all
clc
tEnd = 3         ;%�q0�}�l�p��tEnd�� 
dt   = 0.01     ;%�Cdt��p��@��
x0   = -pi        ;%X��l��m
xEnd = pi         ;%X������m
dx   = 0.1     ;%�Cdx���@��
a    = 0.5       ;%constan(�i�t)
cfl    = a*dt/dx   ;%CFL number
type = 1         ;%��l����A1��step func�A2��sin wave
%%
x = x0 : dx : xEnd;%���Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����
%%
switch type
    case 1  %setp func
        u = heaviside(x);%��l����
        uu = heaviside(x);%�����
    case 2  %sin wave
        u = sin(x);%��l����
        uu   = sin(x);%�����
end
ap = max(a,0);
am = min(a,0);

for n = 2 : length(t)%�����ɶ�(t)����

    u = scheme(u,cfl,1);
	
    %exact part
    switch type
    case 1  %setp func
        uu = heaviside(x-a*dt*n);%�����
        u(1) = 0;%von Newmann BC
    case 2  %sin wave
        
        uu   = sin(x-a*dt*n);%�����
        u(1) = u(end);%periodic BC
end
    
    plot(x,u,'*',x,uu,'--')
    axis([x0, xEnd, min(uu), max(uu)])
    grid on
    xlabel('x');%�����y�ЦW��
    ylabel('u(t)');%�����y�ЦW��
    title(['time(t) = ',num2str(t(n))]); % �ϧμ��D
    pause (dt)
end