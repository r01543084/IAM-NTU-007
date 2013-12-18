%%
%Semi-classical ES-BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007
%%
%��J����
clear all;clc;
x0     = 0            ;% X��l��m
xEnd   = 1            ;% X������m
dx     = 0.1          ;% �Cdx���@��
dt     = 0.1          ;% �Cdt���@��
tEnd   = 30           ;% �q0�}�l�p��tEnd��
relax_time = 1/10000  ;% Relaxation time
cfl    = 1            ;% CFL number

%%
%�ɶ� and ��m����
x = x0:dx:xEnd;
time = 0:dt:tEnd;
nx = length(x);

%�t�תŶ�����(Velocity-Space)
nv = 60;%�]�� Gauss-Hermite ��nv���I�A���F�n���t��domain
        %order �� 2*nv-1
[micro_v,weigt] = GaussHermite(nv);
micro_v = repmat(micro_v,1,nx);%�N���[�t�צV���t�תŶ����i�}
weigt = repmat(weigt,1,nx);%�Nweight���]�V���t�תŶ����i�}�A�H�Q�n��

%%
%initial condition
% Load Fugacity, Macroscopic Velocity and Temperature
% depend on x length
[z0,u0,T0] = IC(nx);

%�N�U���[�q�V���t�תŶ����i�}
z = repmat(z0,nv,1); ux = repmat(u0,nv,1); T = repmat(T0,nv,1);
