clear all
clc
close all
x=-3*pi:0.5:3*pi;
y=-3*pi:0.5:3*pi;

u = 2;
v = 1;
dt=0.05;
t=0:dt:3;
for i = 1:length(t)

[X,Y] = meshgrid(x,y);
U = heaviside(X-u*dt*i)-heaviside(X-u*dt*i-2)+heaviside(Y-v*dt*i)-heaviside(Y-2-v*dt*i);
mesh(x,y,U)
pause(dt)
end