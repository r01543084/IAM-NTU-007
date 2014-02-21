%%
%stability �Ҽ��g��linear advection equation (velocity >0)
%Atmosphere @NTU 2014.1.6
%%
%��J����
clear all
close all
clc
x0   = -2*pi        ;%X��l��m
xEnd = 2*pi         ;%X������m
dx   = 0.1          ;%�Cdx���@��
tEnd = 0.1           ;%�q0�}�l�p��tEnd�� 
cfl  = 0.04          ;%CFL number
a    = 10            ;%constant
dt   = cfl*dx/a     ;%�Cdt���@��
scheme_type = 6    ;%��h���
inital_type = 1    ;%��l����A1��step func�A2��sin wave
%%
x = x0 : dx : xEnd;%���Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����
%%
switch inital_type
    case 1  %setp func
        u  = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%��l����
        uu = exp(-(x+4).^2)+heaviside(x)-heaviside(x-3);%�����
    case 2  %sin wave
        u  = sin(x);%��l����
        uu = sin(x);%�����
end

%%
if scheme_type == 6
    aj = ones(1,length(x)-3)*cfl/4;
    bj = -ones(1,length(x)-3)*cfl/4;
    dj = ones(1,length(x)-2);
else
    aj=0;bj=0;dj=0;%��a b c d ����
end

%% Main Loop
for n = 2 : length(t)%�����ɶ�(t)����
    
    %mean equation
    u = scheme(u,aj,bj,dj,cfl,scheme_type);
	
    %exact part
    switch inital_type
    case 1  %setp func
        uu = exp(-(x+4-a*dt*n).^2)+heaviside(x-a*dt*n)-heaviside(x-3-a*dt*n);%�����
        u(1) = 0;%von Newmann BC
        u(end) = 0;
    case 2  %sin wave
        uu   = sin(x-a*dt*n);%�����
        u(end-5:end) = uu(end-5:end);%periodic BC
        u(1:5) = uu(1:5);
    end
    
    plot(x,u,'-O',x,uu,'--')
    axis([x0, xEnd, min(uu)-0.5, max(uu)+0.5])
    grid on
    xlabel('x');%�����y�ЦW��
    ylabel('u(t)');%�����y�ЦW��
    title(['time(t) = ',num2str(t(n))]); % �ϧμ��D
	pause (dt)
end