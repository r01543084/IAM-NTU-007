%%
%Semi-classical ES-BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007

%% 輸入條件
clear all;close all;clc;
tic
x0     = 0            ;% X初始位置
xEnd   = 1            ;% X結束位置
dx     = 0.01          ;% 每dx切一格
tEnd   = 0.1           ;% 從0開始計算tEnd秒
relax_time = 1/1000  ;% Relaxation time
cfl    = 0.1            ;% CFL number
dv     = 3              ;%polytropic constant
global gamma
gamma  = (dv+2)/dv;           
[xx,dexact,uexact,pexact,mexact,entroexact,energexact] = EulerExact(1,0,1,0.125,0,0.1,tEnd);
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
mirco_v = repmat(mirco_v,1,nx);%將微觀速度向”速度空間”展開
weight  = repmat(weight ,1,nx);%將weight項也向”速度空間”展開，以利積分

%時間離散
dt = dx*cfl/max(mirco_v(:,1))
time = 0:dt:tEnd;

%% 輸入初始條件
%initial condition
% Load Fugacity, Macroscopic Velocity and Temperature
% depend on x length
[marco_u0,T0,density0] = IC(nx);

%將各巨觀量向”速度空間”展開
density = repmat(density0,nv,1);%density
marco_u = repmat(marco_u0,nv,1);%Macroscopic Velocity 
T       = repmat(T0      ,nv,1);%Temperature

%將各巨觀量帶入平衡態方程式
f0 = f_equilibrium(marco_u,mirco_v,T,density,dv);%平衡態方程式

%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)
[density1,marco_u,T] = densityfunc(f0,weight,marco_u,mirco_v,T,density,dv,nv);

%% Marching Scheme
f = f0;% Load initial condition
                 
ap = max(mirco_v,0);%u+
am = min(mirco_v,0);%u-
for tstep = time
    %在每個時間步中，利用各微觀量，更新平衡態分布函數
    f_eq = f_equilibrium(marco_u,mirco_v,T,density,dv);

    %為了與conservation law的符號相等，故使用u
	u_eq = f_eq;
	u = f;
    
     u(:,2:end-1) = u(:,2:end-1) - dt/dx*ap(:,1:end-2).*(u(:,2:end-1)-u(:,1:end-2))...
                                - dt/dx*am(:,2:end-1).*(u(:,3:end)-u(:,2:end-1))...
                                +(dt/relax_time)*(u_eq(:,2:end-1)-u(:,2:end-1));
    %main equation
%     [uR, uL] = weno3(u);
%     u = u - dt/dx*LF_flux(mirco_v,uR,uL);
    f = u;
    
    %利用新得到的f(分布函數)求得可得數密度(n)、通量or動量密度(j_x or nu)、能量密度
    [density,marco_u,T] = densityfunc(f,weight,marco_u,mirco_v,T,density,dv,nv);
    
    % UPDATE macroscopic properties 
	% (here lies a paralellizing computing chalenge)
	%[marco_u,p,energy] = macroproperties1d(density,marco_u,T,nv,dv);
 
    %plot part
    subplot(1,3,1); 
    plot(x,density(1,:),'.',xx,dexact,'--');
    axis tight; title('Density')
    grid on
    subplot(1,3,2)
mesh(f_eq)
view(135,45)
    subplot(1,3,3)
mesh(f)
view(135,45)

%     subplot(2,3,2); plot(x,p(1,:),'.',xx,pexact,'--'); 
%     axis tight; title('Pressure')
%     grid on
%     subplot(2,3,3); plot(x,T(1,:),'.'); 
%     axis tight; title('Temperature')
%     grid on
%     subplot(2,3,5); plot(x,marco_u(1,:),'.',xx,uexact,'--');
%     axis tight; title('velocity in x')
%     grid on
% %     subplot(2,3,6); plot(x,energy(1,:),'.',xx,energexact,'--');
% %     axis tight; title('Energy Density')
% %     grid on
      pause(dt)
end



toc








