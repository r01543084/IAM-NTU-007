%%
%linerAdvection use upwind by Atmosphere at NTU. 2013,12,3
clear all
clc
close all
%%
%initial condition
x0 = -3*pi;
xEnd = 3*pi;
dx = 0.3;
y0 = -3*pi;
yEnd = 3*pi;
dy = 0.3;
tEnd = 3;
u = 10;%x方向速度
v = 10;%y方向速度
cflx = .1;
cfly = .1;
dt = min(cflx*dx/u,cfly*dy/v);
cflx = dt*u/dx
cfly = dt*v/dy

%%
%時間空間離散
x= x0:dx:xEnd;
y= y0:dy:yEnd;
t = 0:dt:tEnd;

%%
[X,Y] = meshgrid(x,y);
%U = heaviside(X)-heaviside(X-2)+heaviside(Y)-heaviside(Y-2);
%U = heaviside(X)-heaviside(X-2)+sin(Y);
U = exp(-*(X).^2-(Y).^2);
for i = 1:length(t)
    U_next = zeros(length(x),length(y));
    U_next(:,2:end) = U(:,2:end)-cflx*(U(:,2:end)-U(:,1:end-1));%x
    U_next(:,1) = U(:,end);%x
    U = U_next;
    U_next(2:end,:) = U(2:end,:)-cfly*(U(2:end,:)-U(1:end-1,:));%y
    U_next(1,:) = U(end,:);%y
    U = U_next;
    mesh(x,y,U)
    pause(dt)
end