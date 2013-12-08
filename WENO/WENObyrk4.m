%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%               WENO3   5th order by RK4                                  %
%                      Atmosphere@ntu 2013.12.2                           %
%                                   thanks Tony                           %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%��J����
clear all;close all;clc;
x0   = -2*pi        ;%X��l��m
xEnd = 2*pi         ;%X������m
dx   = 0.1         ;%�Cdx���@��
tEnd = 3            ;%�q0�}�l�p��tEnd�� 
CFL  = 0.01          ;%CFL number
a    = 10            ;%constant
dt   = CFL*dx/a     ;%�Cdt���@��
%%
x = x0 : dx : xEnd;%���Ŷ�����
t = 0  : dt : tEnd;%���ɶ�����
%%
%mean progream
%u = heaviside(x)-heaviside(x-3);%initial condition
u = sin(x);
for j = 1:length(t)
    F=[];%�b��F
    dF=[];%same
    f=[];%same
    u_weno=[];%same
    u_next=[];%�b��u_next
    apha  = max(abs(a));%�L��flux
    u_weno = weno3new(u);%use weno3 method
    u_weno = [u_weno(1,end) u_weno(1,:);
              u_weno(2,:) u_weno(2,end)];%periodic bc
    f = a*u_weno;%flux term
	F = 0.5*((f(2,:)+f(1,:))-apha*(u_weno(2,:)-u_weno(1,:)));%LF flux
	u_next = u - ( F(2:end)-F(1:end-1) )/dx*dt;%mean equation
    %uu = heaviside(x-a*dt*j)-heaviside(x-3-a*dt*j);
    uu = sin(x-a*dt*j);
    u_next(1)=u(end);%bc
    u = u_next;%renew u
    plot(x,u,'.',x,uu,'--')
    axis([x0-1, xEnd+1, min(uu)-0.5, max(uu)+0.5])
    xlabel('x');%�����y�ЦW��
    ylabel('u(t)');%�����y�ЦW��
    title(['time(t) = ',num2str(t(j))]); % �ϧμ��D
    grid on
    pause(dt)
end