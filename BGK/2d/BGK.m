%%  ---------------2D----------
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
y0     = 0              ;% Y��l��m
yEnd   = 1              ;% Y������m
dx     = 0.01          ;% �Cdx���@��
dy     = 0.01          ;% �Cdy���@��
tEnd   = 0.3            ;% �q0�}�l�p��tEnd��
relax_time = 10^-8      ;% Relaxation time
cfl    = 2              ;% CFL nuber
dv     = 3              ;%polytropic constant
ic_case = 5             ;%IC. case
%% �����ɶ��B�t�תŶ��B��m�Ŷ�
% ��m�Ŷ�����
x = x0:dx:xEnd;
nx = length(x);

y = y0:dy:yEnd;
ny = length(y);

%�t�תŶ�����(Velocity-Space)

%x velocity
nvx = 34;%�]�� Gauss-Hermite ��nv���I�A���F�n���t��domain
        %order �� 2*nv-1
[mirco_vx,weightx] = GaussHermite(nvx);%for integrating range: -inf to inf
weightx = weightx.*exp(mirco_vx.^2);%real weight if not, chack out website
% http://www.efunda.com/math/num_integration/findgausshermite.cfm

%y velocity
nvy = 34;
[mirco_vy,weighty] = GaussHermite(nvy);
weighty = weighty.*exp(mirco_vy.^2);

%�ɶ�����
dt = min(dx*cfl/max(abs(mirco_vx)),dy*cfl/max(abs(mirco_vy)));
time = 0:dt:tEnd;
ap_coef = dt/relax_time;
%% ��J��l����ñN��l�����Jexact sol.��
%initial condition
% depend on x length and y length
fprintf('��J��l����...')
[marco_ux,marco_uy,T,density,p] = IC(ic_case,nx,ny);
fprintf('����\n')
% index �i�}
fprintf('index �i�}...')
    %�Хߥѡ��Ŷ����V���t�תŶ����i�}��index
    id_space = repmat(reshape((1:nx*ny),[ny,nx]),[1,1,nvx,nvy]);%�Ŷ��i�}

    %�Хߥѡ��t�תŶ����V���Ŷ����i�}��index
    idvx = reshape(1:nvx  ,[1,1,nvx,1]);%���s�w�q 1:nvx ������
    idvx = repmat(idvx,[ny,nx+1,1,nvy]);%Microscopic Velocity x ������
                                        %nx+1�O�]�� weno �����W�[ x ��V������
    idvy = reshape(1:nvy  ,[1,1,1,nvy]);%���s�w�q 1:nvy ������
    idvy = repmat(idvy,[ny+1,nx,nvx,1]);%Microscopic Velocity y ������
                                        %ny+1�O�]�� weno �����W�[ y ��V������
fprintf('����\n')

%�N�U���[�q�a�J���źA��{��
fprintf('�N�U���[�q�a�J���źA��{��...')
[g0,h0] = f_equilibrium(marco_ux,mirco_vx,marco_uy,mirco_vy,T,density,id_space...
                        ,idvx(:,1:end-1,:,:),idvy(1:end-1,:,:,:));%���źA��{��
fprintf('����\n')
%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)
[density1,marco_ux1,marco_uy1,T1] = densityfunc(g0,h0,weightx,mirco_vx,weighty,mirco_vy...
                                    ,id_space,idvx(:,1:end-1,:,:),idvy(1:end-1,:,:,:));
fprintf('���լO�_�ϥ�BGK�D�o�����[�q�A�P��l����۵�\n')
disp(sum(sum(density-density1)))
disp(sum(sum(marco_ux-marco_ux1)))
disp(sum(sum(marco_uy-marco_uy1)))
disp(sum(sum(T-T1)))

%% Marching Scheme
% Load initial condition
    disp('Loading distribution equation...')
    g = g0;
    h = h0;
    g_eq = g0;
    h_eq = h0;
    disp('Successful')
    
toc
for tstep = time
    clc
    tstep
    %main loop use RK AP. 4th 
    tic
    %x dir
        %AP. RK step 1
        disp('AP. RK x dir. step 1')
        k1 = ( (-LF_flux(mirco_vx,weno3(g,[0,1]),idvx,nvx))        /dx );
        k2 = ( (-LF_flux(mirco_vx,weno3(h,[0,1]),idvx,nvx))        /dx );

        %AP. RK step 2
        disp('AP. RK x dir. step 2')
        k3 = ( (-LF_flux(mirco_vx,weno3(g+k1/2*dt,[0,1]),idvx,nvx))/dx );
        k4 = ( (-LF_flux(mirco_vx,weno3(h+k2/2*dt,[0,1]),idvx,nvx))/dx );

        %AP. RK step 3
        disp('AP. RK x dir. step 3')
        k5 = ( (-LF_flux(mirco_vx,weno3(g+k3/2*dt,[0,1]),idvx,nvx))/dx );
        k6 = ( (-LF_flux(mirco_vx,weno3(h+k4/2*dt,[0,1]),idvx,nvx))/dx );

        %AP. RK step 4
        disp('AP. RK x dir. step 4')
        k7 = ( (-LF_flux(mirco_vx,weno3(g+k5*dt,[0,1]),idvx,nvx))  /dx );
        k8 = ( (-LF_flux(mirco_vx,weno3(h+k6*dt,[0,1]),idvx,nvx))  /dx );

        %AP. RK final step, sum it!
        g_star = g + 1/6*(k1+2*k3+2*k5+k7)*dt;%mean equation
        h_star = h + 1/6*(k2+2*k4+2*k6+k8)*dt;%mean equation
        toc
    %y dir
        %AP. RK step 1
        disp('AP. RK y dir. step 1')
        k1 = ( (-LF_flux(mirco_vy,weno3(g_star,[1,0]),idvy,nvx))        /dy );
        k2 = ( (-LF_flux(mirco_vy,weno3(h_star,[1,0]),idvy,nvx))        /dy );

        %AP. RK step 2
        disp('AP. RK y dir. step 2')
        k3 = ( (-LF_flux(mirco_vy,weno3(g_star+k1/2*dt,[1,0]),idvy,nvx))/dy );
        k4 = ( (-LF_flux(mirco_vy,weno3(h_star+k2/2*dt,[1,0]),idvy,nvx))/dy );

        %AP. RK step 3
        disp('AP. RK y dir. step 3')
        k5 = ( (-LF_flux(mirco_vy,weno3(g_star+k3/2*dt,[1,0]),idvy,nvx))/dy );
        k6 = ( (-LF_flux(mirco_vy,weno3(h_star+k4/2*dt,[1,0]),idvy,nvx))/dy );

        %AP. RK step 4
        disp('AP. RK y dir. step 4')
        k7 = ( (-LF_flux(mirco_vy,weno3(g_star+k5*dt,[1,0]),idvy,nvx))  /dy );
        k8 = ( (-LF_flux(mirco_vy,weno3(h_star+k6*dt,[1,0]),idvy,nvx))  /dy );

        %AP. RK final step, sum it
        disp('AP. RK final step')
        g_star = g_star + 1/6*(k1+2*k3+2*k5+k7)*dt;%mean equation
        h_star = h_star + 1/6*(k2+2*k4+2*k6+k8)*dt;%mean equation
        
        
        disp('AP. RK Successful')
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
    [density,marco_ux,marco_uy,T] = densityfunc(g_star,h_star,weightx,mirco_vx,weighty,mirco_vy...
                                    ,id_space,idvx(:,1:end-1,:,:),idvy(1:end-1,:,:,:));
    
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
	[g_eq,h_eq] = f_equilibrium(marco_ux,mirco_vx,marco_uy,mirco_vy,T,density,id_space...
                        ,idvx(:,1:end-1,:,:),idvy(1:end-1,:,:,:));
                    
                    
    %main loop use RK 4th

        %RK final step, sum it
        disp('RK final step')
        g = (g_star + ap_coef*(g_eq)) / (1+ap_coef);%mean equation
        h = (h_star + ap_coef*(h_eq)) / (1+ap_coef);%mean equation
        
        
        disp('RK Successful')
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
    [density,marco_ux,marco_uy,T,e,p] = densityfunc(g,h,weightx,mirco_vx,weighty,mirco_vy...
                                    ,id_space,idvx(:,1:end-1,:,:),idvy(1:end-1,:,:,:));
    
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
	[g_eq,h_eq] = f_equilibrium(marco_ux,mirco_vx,marco_uy,mirco_vy,T,density,id_space...
                        ,idvx(:,1:end-1,:,:),idvy(1:end-1,:,:,:));  
                    contour(x,y,density,25)
                    axis equal
                    axis([x0 xEnd y0 yEnd])
                    pause(dt)

end



toc