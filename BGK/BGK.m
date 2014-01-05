%%
%Semi-classical ES-BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007

%% ��J����
clear all;close all;clc;
tic
x0     = 0            ;% X��l��m
xEnd   = 1            ;% X������m
dx     = 0.01          ;% �Cdx���@��
tEnd   = 0.1           ;% �q0�}�l�p��tEnd��
relax_time = 1/1000  ;% Relaxation time
cfl    = 0.1            ;% CFL number
dv     = 3              ;%polytropic constant
global gamma
gamma  = (dv+2)/dv;           
[xx,dexact,uexact,pexact,mexact,entroexact,energexact] = EulerExact(1,0,1,0.125,0,0.1,tEnd);
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
mirco_v = repmat(mirco_v,1,nx);%�N�L�[�t�צV���t�תŶ����i�}
weight  = repmat(weight ,1,nx);%�Nweight���]�V���t�תŶ����i�}�A�H�Q�n��

%�ɶ�����
dt = dx*cfl/max(mirco_v(:,1))
time = 0:dt:tEnd;

%% ��J��l����
%initial condition
% Load Fugacity, Macroscopic Velocity and Temperature
% depend on x length
[marco_u0,T0,density0] = IC(nx);

%�N�U���[�q�V���t�תŶ����i�}
density = repmat(density0,nv,1);%density
marco_u = repmat(marco_u0,nv,1);%Macroscopic Velocity 
T       = repmat(T0      ,nv,1);%Temperature

%�N�U���[�q�a�J���źA��{��
f0 = f_equilibrium(marco_u,mirco_v,T,density,dv);%���źA��{��

%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)
[density1,marco_u,T] = densityfunc(f0,weight,marco_u,mirco_v,T,density,dv,nv);

%% Marching Scheme
f = f0;% Load initial condition
                 
ap = max(mirco_v,0);%u+
am = min(mirco_v,0);%u-
for tstep = time
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
    f_eq = f_equilibrium(marco_u,mirco_v,T,density,dv);

    %���F�Pconservation law���Ÿ��۵��A�G�ϥ�u
	u_eq = f_eq;
	u = f;
    
     u(:,2:end-1) = u(:,2:end-1) - dt/dx*ap(:,1:end-2).*(u(:,2:end-1)-u(:,1:end-2))...
                                - dt/dx*am(:,2:end-1).*(u(:,3:end)-u(:,2:end-1))...
                                +(dt/relax_time)*(u_eq(:,2:end-1)-u(:,2:end-1));
    %main equation
%     [uR, uL] = weno3(u);
%     u = u - dt/dx*LF_flux(mirco_v,uR,uL);
    f = u;
    
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
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








