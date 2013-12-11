%%
%Hopscotch method solve 2D heat equation
%by Atmosphere @ NTU 2013.12.11
%%
clear all;close all;clc;
%丁丁瞒床
dx = 0.5;
dy = 0.5;
dt = 0.001;
x0 = 0;
xEnd = 10;
y0 = 0;
yEnd = 10;
tEnd = 10;
k = 50;%heat coef.
x = x0:dx:xEnd;
y = y0:dy:yEnd;
t = 0:dt:tEnd;
u = zeros(length(x),length(y));%initial condition
rx = dt*k/dx^2
ry = dt*k/dy^2
%bc
u(1,:)   = 0;%y=0放
u(end,:) = 20;%y=yEnd放

u(:,1)   = 0;%x=0放
u(:,end) = 10;%x=xEnd放

%main loop
for n = 1:length(t)
    for i = 2:length(x)-1
        for j = 2:length(y)-1
            if mod(i+j+n,2) == 0
                u(i,j) = u(i,j)+(u(i+1,j)-2*u(i,j)+u(i-1,j))*rx ...
                               +(u(i,j+1)-2*u(i,j)+u(i,j-1))*ry;
            else
                u(i,j) = u(i,j)+(u(i+1,j)-2*u(i,j)+u(i-1,j))*rx ...
                               +(u(i,j+1)-2*u(i,j)+u(i,j-1))*ry;
            end
        end
    end
	[C h]=contourf(u);
    clabel(C,h)
    caxis([0 20])
    contourcbar
    pause(dt)
    
end