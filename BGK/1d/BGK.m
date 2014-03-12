%%
%BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007
% thanks Tony dimensionless coef.
%% ��J����
clear all;close all;clc;
tic
x0     = 0              ;% X��l��m
xEnd   = 1              ;% X������m
dx     = 0.005          ;% �Cdx���@��
tEnd   = 0.1            ;% �q0�}�l�p��tEnd��
relax_time = 0.001      ;% Relaxation time
cfl    = 1              ;% CFL nuber
dv     = 3              ;% polytropic constant
apha   = 1              ;% �Ϯg(0)>>>>�l��(1)
global gamma
gamma  = (dv+2)/dv;        

%% �����ɶ��B�t�תŶ��B��m�Ŷ�
% ��m�Ŷ�����
x = x0:dx:xEnd;
nx = length(x);

%�t�תŶ�����(Velocity-Space)
nv = 60;%�]�� Gauss-Hermite ��nv���I�A���F�n���t��domain
        %order �� 2*nv-1
[mirco_v,weight] = GaussHermite(nv);%for integrating range: -inf to inf
weight = weight.*exp(mirco_v.^2);%real weight if not, chack out website
% http://www.efunda.com/math/num_integration/findgausshermite.cfm
weightb = weight;%for boundary weight
weightb(1:ceil(nv/2))=0;%for boundary weight
mirco_vb = (mirco_v+abs(mirco_v))/2;%for boundary mirco velocity
%�ɶ�����
dt = dx*cfl/max(abs(mirco_v));%limit is dt = relax_time*0.8
time = 0:dt:tEnd;

%% ��J��l����ñN��l�����Jexact sol.��
%initial condition
% depend on x length
[marco_u,T,density,p] = IC(nx);

%���ѳ���
[xx,dexact,uexact,pexact,mexact,entroexact,energexact,Texact] = ...
    EulerExact(density(1),marco_u(1),p(1),...
               density(end),marco_u(end),p(end),tEnd);

%�N�U���[�q�V���t�תŶ����i�}
idx = repmat(1:nx,nv,1);%density
idv = repmat((1:nv)',1,nx);%Macroscopic Velocity 

%�N�U���[�q�a�J���źA��{��
[g0,h0] = f_equilibrium(marco_u,mirco_v,T,density,idx,idv);%���źA��{��

%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)
[density1,marco_u1,T1] = densityfunc(g0,h0,weight,mirco_v,idv);
fprintf('���լO�_�ϥ�BGK�D�o�����[�q�A�P��l����۵�\n')
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
    %RK step 1
    k1 = ( (-LF_flux(mirco_v,weno3(g),idv,nv))        /dx );
    k2 = ( (-LF_flux(mirco_v,weno3(h),idv,nv))        /dx );
    
    %RK step 2
    k3 = ( (-LF_flux(mirco_v,weno3(g+k1/2*dt),idv,nv))/dx );
    k4 = ( (-LF_flux(mirco_v,weno3(h+k2/2*dt),idv,nv))/dx );
    
    %RK step 3
    k5 = ( (-LF_flux(mirco_v,weno3(g+k3/2*dt),idv,nv))/dx );
    k6 = ( (-LF_flux(mirco_v,weno3(h+k4/2*dt),idv,nv))/dx );
    
    %RK step 4
    k7 = ( (-LF_flux(mirco_v,weno3(g+k5*dt),idv,nv))  /dx );
    k8 = ( (-LF_flux(mirco_v,weno3(h+k6*dt),idv,nv))  /dx );
    
    %RK final step, sum it!
    g = g + 1/6*(k1+2*k3+2*k5+k7)*dt+(dt/relax_time)*(g_eq-g);%mean equation
    h = h + 1/6*(k2+2*k4+2*k6+k8)*dt+(dt/relax_time)*(h_eq-h);%mean equation
	
    %reflaction and diffusive boundary
    [g(:,end-1),h(:,end-1),g(:,end),h(:,end)] = ...
        boundaryf(g(:,end-1),h(:,end-1),mirco_v,weightb...
                                        ,nv,apha);
    
    
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
    [density,marco_u,T,E,p] = densityfunc(g,h,weight,mirco_v,idv);
   
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
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
    axis tight; title('Density')
    grid on
    subplot(2,3,3); plot(x,p,'.',xx,pexact,'--'); 
    axis tight; title('Pressure')
    grid on
    subplot(2,3,5); plot(x,T,'.',xx,Texact,'--'); 
    axis tight; title('Temperature')
    grid on
    subplot(2,3,6); plot(x,marco_u,'.',xx,uexact,'--');
    axis tight; title('velocity in x')
    grid on
	pause(dt)

end

toc


