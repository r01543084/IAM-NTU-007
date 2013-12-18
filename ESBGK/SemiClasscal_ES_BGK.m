%%
%Semi-classical ES-BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007
%%
%輸入條件
clear all;clc;
x0     = 0            ;% X初始位置
xEnd   = 1            ;% X結束位置
dx     = 0.1          ;% 每dx切一格
dt     = 0.1          ;% 每dt切一格
tEnd   = 30           ;% 從0開始計算tEnd秒
relax_time = 1/10000  ;% Relaxation time
cfl    = 1            ;% CFL number

%%
%時間 and 位置離散
x = x0:dx:xEnd;
time = 0:dt:tEnd;
nx = length(x);

%速度空間離散(Velocity-Space)
nv = 60;%因為 Gauss-Hermite 取nv個點，為了積分速度domain
        %order 為 2*nv-1
[micro_v,weigt] = GaussHermite(nv);
micro_v = repmat(micro_v,1,nx);%將巨觀速度向”速度空間”展開
weigt = repmat(weigt,1,nx);%將weight項也向”速度空間”展開，以利積分

%%
%initial condition
% Load Fugacity, Macroscopic Velocity and Temperature
% depend on x length
[z0,u0,T0] = IC(nx);

%將各巨觀量向”速度空間”展開
z = repmat(z0,nv,1); ux = repmat(u0,nv,1); T = repmat(T0,nv,1);
