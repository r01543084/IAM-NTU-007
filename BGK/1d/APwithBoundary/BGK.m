%%
%BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007
% thanks Tony dimensionless coef.
%% 輸入條件
clear all;close all;clc;
tic
x0     = 0            ;% X初始位置
xEnd   = 1            ;% X結束位置
dx     = 0.01          ;% 每dx切一格
tEnd   = 0.55           ;% 從0開始計算tEnd秒
relax_time = 10^-10  ;% Relaxation time
cfl    = 1            ;% CFL nuber
dv     = 3              ;%polytropic constant
tmethod = 'SSPRK4'      ;%時間方法 SSPRK4 or RK4
apha   = 0              ;% reflaction(0)>>>>diffusitive(1)
global gamma
gamma  = (dv+2)/dv;           

%% 離散時間、速度空間、位置空間
% 位置空間離散
x = x0:dx:xEnd;
nx = length(x);

%速度空間離散(Velocity-Space)
nv = 100;%因為 Gauss-Hermite 取nv個點，為了積分速度domain
        %order 為 2*nv-1
[mirco_v,weight] = GaussHermite(nv);%for integrating range: -inf to inf
weight = weight.*exp(mirco_v.^2);%real weight if not, chack out website
% http://www.efunda.com/math/num_integration/findgausshermite.cfm
weightb = weight;%for boundary weight
weightb(1:ceil(nv/2))=0;%for boundary weight
mirco_vb = (mirco_v+abs(mirco_v))/2;%for boundary mirco velocity

%時間離散
dt = dx*cfl/max(abs(mirco_v));%limit is dt = relax_time*0.8
time = 0:dt:tEnd;
ap_coef = (relax_time+dt)/relax_time;

%% 輸入初始條件並將初始條件輸入exact sol.中
%initial condition
% depend on x length
[marco_u,T,density,p] = IC(nx);

%正解部分
[xx,dexact,uexact,pexact,mexact,entroexact,energexact,Texact] = ...
    EulerExact(density(1),marco_u(1),p(1),...
               density(end),marco_u(end),p(end),tEnd);

%將各巨觀量向”速度空間”展開
idx = repmat(1:nx,[nv,1]);%density
idv = repmat((1:nv)',[1,nx]);%Velocity x domain 

%將各巨觀量帶入平衡態方程式
[g0,h0] = f_equilibrium(marco_u,mirco_v,T,density,idx,idv);%平衡態方程式

%利用平衡態方程式積分(Gauss-Hermite)可得數密度(n)、通量or動量密度(j_x or nu)、
%                                   能量密度(Energy Density)
[density1,marco_u1,T1] = densityfunc(g0,h0,weight,mirco_v,idv);
fprintf('測試是否使用BGK求得之巨觀量，與初始條件相等\n')
disp(sum(density-density1))
disp(sum(marco_u-marco_u1))
disp(sum(T-T1))

%% Marching Scheme
% Load initial condition
g = g0;
h = h0;
g_eq = g0;
h_eq = h0;

for tstep = time
    %main loop use RK 4th
%--AP start

	%for near boundary node
    g_b = dfdx(g(:,end-5:end),mirco_v,idv(:,end-1:end),nv,dx,dt,tmethod,'b');
    h_b = dfdx(h(:,end-5:end),mirco_v,idv(:,end-1:end),nv,dx,dt,tmethod,'b');
    
    %RK final step, sum it!
    g1 = dfdx(g,mirco_v,idv,nv,dx,dt,tmethod,'n');%mean equation
    h1 = dfdx(h,mirco_v,idv,nv,dx,dt,tmethod,'n');%mean equation
    
	%put boundary node info. into near wall 
	g(:,end-1:end) = g_b;
    h(:,end-1:end) = h_b;

    %利用新得到的f(分布函數)求得可得數密度(n)、通量or動量密度(j_x or nu)、能量密度
    [density1,marco_u1,T1] = densityfunc(g1,h1,weight,mirco_v,idv);
    
    %在每個時間步中，利用各微觀量，更新平衡態分布函數
	[g_eq,h_eq] = f_equilibrium(marco_u1,mirco_v,T1,density1,idx,idv);

%--AP end    
    
    
    %RK final step, sum it!
    g = (g1 +(dt/relax_time)*(g_eq)) / ap_coef;%mean equation
    h = (h1 +(dt/relax_time)*(h_eq)) / ap_coef;%mean equation

    %reflaction and diffusive boundary
    [g(:,end-1),h(:,end-1),g(:,end),h(:,end)] = ...
        boundaryf(g(:,end-1),h(:,end-1),mirco_v,weightb,nv,apha);

    %利用新得到的f(分布函數)求得可得數密度(n)、通量or動量密度(j_x or nu)、能量密度
    [density,marco_u,T,e,p] = densityfunc(g,h,weight,mirco_v,idv);
    
    %在每個時間步中，利用各微觀量，更新平衡態分布函數
	[g_eq,h_eq] = f_equilibrium(marco_u,mirco_v,T,density,idx,idv);
        
	%plot part
    subplot(1,3,1); 
    meshc(x,mirco_v,g);
    axis tight; title('bgk function');
    xlabel('x - Spatial Domain');
	ylabel('v - Velocity Space');
	zlabel('f - Probability');
    grid on
    view(135,45)
    subplot(2,3,2); 
    plot(x,density,'.',xx,dexact,'--');
    axis([0,1,0,1.5]); title('Density')
    grid on
    subplot(2,3,3); plot(x,p,'.',xx,pexact,'--'); 
    axis([0,1,0,1]); title('Pressure')
    grid on
    subplot(2,3,5); plot(x,T,'.',xx,Texact,'--'); 
    axis([0,1,0.5,2.3]);title('Temperature')
    grid on
    subplot(2,3,6); plot(x,marco_u,'.',xx,uexact,'--');
    axis([0,1,-0.1,1]); title('velocity in x')
    grid on
	pause(dt)   

end



toc
