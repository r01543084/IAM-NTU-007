%%
%Semi-classical ES-BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007

%% 輸入條件
clear all;clc;
pause()
x0     = 0            ;% X初始位置
xEnd   = 1            ;% X結束位置
dx     = 0.01          ;% 每dx切一格
tEnd   = 0.1           ;% 從0開始計算tEnd秒
relax_time = 1/10000  ;% Relaxation time
cfl    = 0.1            ;% CFL number
theta  = -1           ;% (-1) BE, (0) MB, (1) FD.
                       
%% 離散時間、速度空間、位置空間
% 位置空間離散
x = x0:dx:xEnd;
nx = length(x);
 
%速度空間離散(Velocity-Space)
nv = 60;%因為 Gauss-Hermite 取nv個點，為了積分速度domain
        %order 為 2*nv-1
[mirco_v,weight] = GaussHermite(nv);%for integrating range: -inf to inf
weight = weight.*exp(mirco_v.^2);%real weight if not, chack out website
% http://www.efunda.com/math/num_integration/findgausshermite.cfm
mirco_v = repmat(mirco_v,1,nx);%將巨觀速度向”速度空間”展開
weight  = repmat(weight ,1,nx);%將weight項也向”速度空間”展開，以利積分

%時間離散
dt = dx*cfl/max(mirco_v(:,1));
time = 0:dt:tEnd;

%% 輸入初始條件
%initial condition
% Load Fugacity, Macroscopic Velocity and Temperature
% depend on x length
[z0,marco_u0,T0] = IC(nx);

%將各巨觀量向”速度空間”展開
z       = repmat(z0      ,nv,1);%Fugacity
marco_u = repmat(marco_u0,nv,1);%Macroscopic Velocity 
T       = repmat(T0      ,nv,1);%Temperature

%將各巨觀量帶入平衡態方程式
f0 = f_equilibrium(z,marco_u,mirco_v,T,theta);%平衡態方程式

%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)
[n,j_x,epsilon] = densityfunc(weight,f0,mirco_v);

%% Marching Scheme
f = f0;% Load initial condition

for tstep = time
    %在每個時間步中，利用各微觀量，更新平衡態分布函數
    f_eq = f_equilibrium(z,marco_u,mirco_v,T,theta);
    
    %為了與conservation law的符號相等，故使用u
	u_eq = f_eq;
	u = f;
    
    %main equation
    [uR, uL] = weno3(u);
    u = u - dt/dx*LF_flux(mirco_v,uR,uL)+(dt/relax_time)*(u_eq-u);
    f = u;
    
    %利用新得到的f(分布函數)求得可得數密度(n)、通量or動量密度(j_x or nu)、能量密度
    [n,j_x,epsilon] = densityfunc(weight,f,mirco_v);
    
    % UPDATE macroscopic properties 
	% (here lies a paralellizing computing chalenge)
	[z,marco_u,T,pressure] = macroproperties1d(n,j_x,epsilon,nx,nv,theta);
    
    %plot part
    subplot(2,3,1); plot(x,n(1,:),'.');
    axis tight; title('Density')
    grid on
    subplot(2,3,2); plot(x,pressure(1,:),'.'); 
    axis tight; title('Pressure')
    grid on
    subplot(2,3,3); plot(x,T(1,:),'.'); 
    axis tight; title('Temperature')
    grid on
    subplot(2,3,4); plot(x,z(1,:),'.'); 
    axis tight; title('Fugacity')
    grid on
    subplot(2,3,5); plot(x,marco_u(1,:),'.');
    axis tight; title('velocity in x')
    grid on
    subplot(2,3,6); plot(x,epsilon(1,:),'.');
    axis tight; title('Energy Density')
    grid on
    pause(dt)
end














