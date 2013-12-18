%%
%Semi-classical ES-BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007

%% ��J����
clear all;clc;
x0     = 0            ;% X��l��m
xEnd   = 1            ;% X������m
dx     = 0.1          ;% �Cdx���@��
tEnd   = 30           ;% �q0�}�l�p��tEnd��
relax_time = 1/10000  ;% Relaxation time
cfl    = 1            ;% CFL number
theta  = -1           ;% (-1) BE, (0) MB, (1) FD.
                       
%% �����ɶ��B�t�תŶ��B��m�Ŷ�
% ��m�Ŷ�����
x = x0:dx:xEnd;
nx = length(x);

%�t�תŶ�����(Velocity-Space)
nv = 60;%�]�� Gauss-Hermite ��nv���I�A���F�n���t��domain
        %order �� 2*nv-1
[micro_v,weight] = GaussHermite(nv);%for integrating range: -inf to inf
weight = weight.*exp(micro_v.^2);%real weight if not, chack out website
% http://www.efunda.com/math/num_integration/findgausshermite.cfm
micro_v = repmat(micro_v,1,nx);%�N���[�t�צV���t�תŶ����i�}
weight  = repmat(weight ,1,nx);%�Nweight���]�V���t�תŶ����i�}�A�H�Q�n��

%�ɶ�����
dt = dx*CFL/max(micro_v(:,1));
time = 0:dt:tEnd;

%% ��J��l����
%initial condition
% Load Fugacity, Macroscopic Velocity and Temperature
% depend on x length
[z0,marco_u0,T0] = IC(nx);

%�N�U���[�q�V���t�תŶ����i�}
z       = repmat(z0      ,nv,1);%Fugacity
marco_u = repmat(marco_u0,nv,1);%Macroscopic Velocity 
T       = repmat(T0      ,nv,1);%Temperature

%�N�U���[�q�a�J���źA��{��
f0 = f_equilibrium(z,marco_u,micro_v,T,theta);%���źA��{��

%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)
[n,j_x,epsilon] = densityfunc(weight,f0,micro_v);

%% Marching Scheme
f = f0;% Load initial condition
for tstep = time
    
    
    
    
    
    %plot part
    subplot(2,3,1); plot(x,n(1,:),'o');
    axis tight; title('Density')
    grid on
    subplot(2,3,2); plot(x,pressure(1,:),'o'); 
    axis tight; title('Pressure')
    grid on
    subplot(2,3,3); plot(x,T(1,:),'o'); 
    axis tight; title('Temperature')
    grid on
    subplot(2,3,4); plot(x,z(1,:),'o'); 
    axis tight; title('Fugacity')
    grid on
    subplot(2,3,5); plot(x,marco_u(1,:),'o');
    axis tight; title('velocity in x')
    grid on
    subplot(2,3,6); plot(x,epsilon(1,:),'o');
    axis tight; title('Energy Density')
    grid on
end













