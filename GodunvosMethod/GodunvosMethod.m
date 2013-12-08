%%
%stability �Ҽ��g��linear advection equationGodunov��s method
%�Ŷ������ϥ�upwind method
%�ɶ������ϥ�euler forward method
%%
%��J����
clear all
clc
tEnd = 1       ;%�q0�}�l�p��10�� 
dt   = 0.0001    ;%�C0.01��p��@��
x0   = -1       ;%X��l��m
xEnd =  5      ;%X������m
dx   = 0.001    ;%�C0.1���@��
%%
x = x0 : dx : xEnd;%���Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����
u = exp(-16*x.^2)+1;%��l����
uCompare = exp(-16*x.^2)+1;%�����
%%
%�D�{��
for j = 1:length(t)
    u_next = zeros(size(u));%�ŧi�s�ܶq�H�K�^��X��
    for i = 2 : length(x)-1
        u_next(i) = u(i) - dt/dx * (F(u(i),u(i+1))-F(u(i-1),u(i)));%mean equation
    end
    u = u_next;%update info
    u(1) = 1;%BC
    plot(x,u,'*',x,uCompare,'--')
    axis([x0, xEnd, min(uCompare), max(uCompare)])
    grid on
    xlabel('x');%�����y�ЦW��
    ylabel('u(t)');%�����y�ЦW��
    title(['time(t) = ',num2str(t(j))]); % �ϧμ��D
    pause(dt)
end