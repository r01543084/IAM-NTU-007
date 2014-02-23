%%  ---------------2D----------
%Classical BGK equation
%Using WENO3 method 
%(fifth order in space domain fourth order in time domain)
%and Gauss-Hermite intergel.
%by Atmosphere @ NTU 2013.12.17
% reference NTU IAM. lab 007
% thanks Tony dimensionless coef.
%% ��J����
clear all;close all;clc;
time_ic = [0.2 0.2  0.3 0.25 0.23 0.3 0.25 0.25 0.3 0.15 ... 
           0.3 0.25 0.3 0.1  0.2  0.2 0.3  0.2  0.3      ];
for ic_case = 7
tic
name      ='CBGK2d';  % Simulation Name
x0     = 0                  ;% X��l��m
xEnd   = 1                  ;% X������m
y0     = 0                  ;% Y��l��m
yEnd   = 1                  ;% Y������m
dx     = 1/200               ;% �Cdx���@��
dy     = 1/200               ;% �Cdy���@��
tEnd   = time_ic(ic_case)   ;% �q0�}�l�p��tEnd��
r_time = 10^-8              ;% Relaxation time
cfl    = 1                  ;% CFL nuber
%ic_case = 7                ;%IC. case
write_ans = 1               ;%(0)no,(1)yes
draw = 2                    ;%(0)no,(1)in time, (2)last time
using_time = 0              ;%how many time we use in this case
pic_num = 20                ;%how many picture u take
%% �����ɶ��B�t�תŶ��B��m�Ŷ�
% ��m�Ŷ�����
x = x0:dx:xEnd;
nx = length(x);

y = y0:dy:yEnd;
ny = length(y);

%�t�תŶ�����(Velocity-Space)

    %x velocity
    nvx = 36;%�]�� Gauss-Hermite ��nv���I�A���F�n���t��domain
            %order �� 2*nv-1
    [mirco_vx,weightx] = GaussHermite(nvx);%for integrating range: -inf to inf
    weightx = weightx.*exp(mirco_vx.^2);%real weight if not, chack out website
    % http://www.efunda.com/math/num_integration/findgausshermite.cfm

    %y velocity
    nvy = 36;
    [mirco_vy,weighty] = GaussHermite(nvy);
    weighty = weighty.*exp(mirco_vy.^2);

%�ɶ�����
    dt = min(dx*cfl/max(abs(mirco_vx)),dy*cfl/max(abs(mirco_vy)));
    time = [0:dt:tEnd-dt tEnd];%�j���̫�@�B��tnd
    ap_coef = dt/r_time;
    shutter_time = ceil((length(time)-1)/pic_num);%�C�h�֨B��@�i�Ӥ�

%% �s�ɮ׳���
    if write_ans == 1
        %�}�ɦW
        [ID, IDn] = ID_name(name,0,nx,ny,nvx,nvy,4,r_time,ic_case);
        %�и�Ƨ�
        mkdir(ID,'T');mkdir(ID,'density');mkdir(ID,'p');
        mkdir(ID,'e');mkdir(ID,'Ux');mkdir(ID,'Uy');
        %save parameter
        all_parameter = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/all_parameter.mat'];
        save(all_parameter);
    end

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
    parfor i = 1:nx
    [g0(:,i,:,:),h0(:,i,:,:)] = f_equilibrium(marco_ux,mirco_vx,marco_uy,mirco_vy,T,density,id_space(:,i,:,:)...
                            ,idvx(:,i,:,:),idvy(1:end-1,i,:,:));%���źA��{��
    end
                    
    fprintf('����\n')
%�Q�Υ��źA��{���n��(Gauss-Hermite)�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B
%                                   ��q�K��(Energy Density)
    parfor i=1:nx
        [density1(:,i),marco_ux1(:,i),marco_uy1(:,i),T1(:,i)]...
            = densityfunc(g0(:,i,:,:),h0(:,i,:,:),weightx,mirco_vx...
            ,weighty,mirco_vy,idvx(:,i,:,:),idvy(1:end-1,i,:,:));
    end
    
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

    %main loop
    counter = 0             ;%counter
for tstep = time
    clc
    tstep
    %main loop use RK AP. 4th
    
    %x dir
    disp('AP. RK x dir. ')
    parfor i = 1:ny
        %AP. RK step 1
        kx1 = ( (-LF_flux(mirco_vx,weno3(g(i,:,:,:),[0,1]),idvx(i,:,:,:),nvx))        /dx );
        kx2 = ( (-LF_flux(mirco_vx,weno3(h(i,:,:,:),[0,1]),idvx(i,:,:,:),nvx))        /dx );

        %AP. RK step 2
        kx3 = ( (-LF_flux(mirco_vx,weno3(g(i,:,:,:)+kx1/2*dt,[0,1]),idvx(i,:,:,:),nvx))/dx );
        kx4 = ( (-LF_flux(mirco_vx,weno3(h(i,:,:,:)+kx2/2*dt,[0,1]),idvx(i,:,:,:),nvx))/dx );

        %AP. RK step 3
        kx5 = ( (-LF_flux(mirco_vx,weno3(g(i,:,:,:)+kx3/2*dt,[0,1]),idvx(i,:,:,:),nvx))/dx );
        kx6 = ( (-LF_flux(mirco_vx,weno3(h(i,:,:,:)+kx4/2*dt,[0,1]),idvx(i,:,:,:),nvx))/dx );

        %AP. RK step 4
        kx7 = ( (-LF_flux(mirco_vx,weno3(g(i,:,:,:)+kx5*dt,[0,1]),idvx(i,:,:,:),nvx))  /dx );
        kx8 = ( (-LF_flux(mirco_vx,weno3(h(i,:,:,:)+kx6*dt,[0,1]),idvx(i,:,:,:),nvx))  /dx );

        %AP. RK final step, sum it!
        g_star(i,:,:,:) = g(i,:,:,:) + 1/6*(kx1+2*kx3+2*kx5+kx7)*dt;%mean equation
        h_star(i,:,:,:) = h(i,:,:,:) + 1/6*(kx2+2*kx4+2*kx6+kx8)*dt;%mean equation
    end
    
    %y dir
    disp('AP. RK y dir.')
    parfor j = 1:nx
        %AP. RK step 1
        ky1 = ( (-LF_flux(mirco_vy,weno3(g(:,j,:,:),[1,0]),idvy(:,j,:,:),nvy))        /dy );
        ky2 = ( (-LF_flux(mirco_vy,weno3(h(:,j,:,:),[1,0]),idvy(:,j,:,:),nvy))        /dy );

        %AP. RK step 2
        ky3 = ( (-LF_flux(mirco_vy,weno3(g(:,j,:,:)+ky1/2*dt,[1,0]),idvy(:,j,:,:),nvy))/dy );
        ky4 = ( (-LF_flux(mirco_vy,weno3(h(:,j,:,:)+ky2/2*dt,[1,0]),idvy(:,j,:,:),nvy))/dy );

        %AP. RK step 3
        ky5 = ( (-LF_flux(mirco_vy,weno3(g(:,j,:,:)+ky3/2*dt,[1,0]),idvy(:,j,:,:),nvy))/dy );
        ky6 = ( (-LF_flux(mirco_vy,weno3(h(:,j,:,:)+ky4/2*dt,[1,0]),idvy(:,j,:,:),nvy))/dy );

        %AP. RK step 4
        ky7 = ( (-LF_flux(mirco_vy,weno3(g(:,j,:,:)+ky5*dt,[1,0]),idvy(:,j,:,:),nvy))  /dy );
        ky8 = ( (-LF_flux(mirco_vy,weno3(h(:,j,:,:)+ky6*dt,[1,0]),idvy(:,j,:,:),nvy))  /dy );

    
        %AP. RK final step, sum it
        g_star(:,j,:,:) = g_star(:,j,:,:) + 1/6*(ky1+2*ky3+2*ky5+ky7)*dt;%mean equation
        h_star(:,j,:,:) = h_star(:,j,:,:) + 1/6*(ky2+2*ky4+2*ky6+ky8)*dt;%mean equation
    end
        
    disp('AP. RK Successful')
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
    parfor i=1:nx
        [density(:,i),marco_ux(:,i),marco_uy(:,i),T(:,i),e(:,i),p(:,i)]...
            = densityfunc(g_star(:,i,:,:),h_star(:,i,:,:),weightx,mirco_vx...
            ,weighty,mirco_vy,idvx(:,i,:,:),idvy(1:end-1,i,:,:));
    end
    
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
    parfor i = 1:nx
        [g_eq(:,i,:,:),h_eq(:,i,:,:)] = f_equilibrium(marco_ux...
            ,mirco_vx,marco_uy,mirco_vy,T,density,id_space(:,i,:,:)...
            ,idvx(:,i,:,:),idvy(1:end-1,i,:,:));%���źA��{��
    end
    
    %main loop use RK 4th
    
    %RK final step, sum it
    disp('RK final step')
    g = (g_star + ap_coef*(g_eq)) / (1+ap_coef);%mean equation
    h = (h_star + ap_coef*(h_eq)) / (1+ap_coef);%mean equation
    
    
    disp('RK Successful')
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
    parfor i=1:nx
        [density(:,i),marco_ux(:,i),marco_uy(:,i),T(:,i),e(:,i),p(:,i)]...
            = densityfunc(g(:,i,:,:),h(:,i,:,:),weightx,mirco_vx...
            ,weighty,mirco_vy,idvx(:,i,:,:),idvy(1:end-1,i,:,:));
    end
    
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
    
    parfor i = 1:nx
        [g_eq(:,i,:,:),h_eq(:,i,:,:)] = f_equilibrium(marco_ux...
            ,mirco_vx,marco_uy,mirco_vy,T,density,id_space(:,i,:,:)...
            ,idvx(:,i,:,:),idvy(1:end-1,i,:,:));%���źA��{��
    end
    
    %Saving part
    if write_ans == 1 && (mod(tstep,shutter_time*dt) == 0 || tstep == time(end))
        counter = counter+1;
        tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/T/T',num2str(counter),'.mat'];
        save(tt,'T','tstep');
        tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/density/density',num2str(counter),'.mat'];
        save(tt,'density','tstep');
        tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/p/p',num2str(counter),'.mat'];
        save(tt,'p','tstep');
        tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/e/e',num2str(counter),'.mat'];
        save(tt,'e','tstep');
        tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/Ux/Ux',num2str(counter),'.mat'];
        save(tt,'marco_ux','tstep');
        tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/Uy/Uy',num2str(counter),'.mat'];
        save(tt,'marco_uy','tstep');
    end
    
    if draw == 1
        %plot
        subplot(2,3,1); contourf(x,y,T,25); axis([x0,xEnd,y0,yEnd]);
        axis equal; title({['Temperature'];['time=',num2str(tstep)]});grid on
        subplot(2,3,2); contourf(x,y,density,25); axis([x0,xEnd,y0,yEnd]);
        axis equal; title({['Density'];['time=',num2str(tstep)]});grid on
        subplot(2,3,3); contourf(x,y,e,25); axis([x0,xEnd,y0,yEnd]);
        axis equal; title({['Internal Energy'];['time=',num2str(tstep)]});grid on
        subplot(2,3,4); contourf(x,y,p,25); axis([x0,xEnd,y0,yEnd]);
        axis equal; title({['Pressure'];['time=',num2str(tstep)]});grid on
        subplot(2,3,5); contourf(x,y,marco_ux,25); axis([x0,xEnd,y0,yEnd]);
        axis equal; title({['Marcoscopic Ux'];['time=',num2str(tstep)]});grid on
        subplot(2,3,6); contourf(x,y,marco_uy,25); axis([x0,xEnd,y0,yEnd]);
        axis equal; title({['Marcoscopic Uy'];['time=',num2str(tstep)]});grid on
        pause(0.5)
    end
end
    
%saving this case using time
using_time = toc;
tt = ['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/using_time','.mat'];
save(tt,'using_time','counter');

%% plot part
    % if you don't have any initial condition, then double click
    % all_parameter. And if no drawing, recheck 'draw' and 'counter' number.
    % And change ur '/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/' to your
    % flie add..
    if draw == 2
        for i =1:counter
            %tell matlab go where to find data
            TT=['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/T/T',num2str(i),'.mat'];
            DD=['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/density/density',num2str(i),'.mat'];
            pp=['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/p/p',num2str(i),'.mat'];
            ee=['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/e/e',num2str(i),'.mat'];
            UUx=['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/Ux/Ux',num2str(i),'.mat'];
            UUy=['/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/',ID,'/Uy/Uy',num2str(i),'.mat'];
            
            %load data
            load(TT);load(DD);load(pp);load(ee);load(UUx);load(UUy);
            
            %plot
            subplot(2,3,1); contourf(x,y,T,25); axis([x0,xEnd,y0,yEnd]);
            axis equal; title({['Temperature'];['time=',num2str(tstep)]});grid on
            subplot(2,3,2); contourf(x,y,density,25); axis([x0,xEnd,y0,yEnd]);
            axis equal; title({['Density'];['time=',num2str(tstep)]});grid on
            subplot(2,3,3); contourf(x,y,e,25); axis([x0,xEnd,y0,yEnd]);
            axis equal; title({['Internal Energy'];['time=',num2str(tstep)]});grid on
            subplot(2,3,4); contourf(x,y,p,25); axis([x0,xEnd,y0,yEnd]);
            axis equal; title({['Pressure'];['time=',num2str(tstep)]});grid on
            subplot(2,3,5); contourf(x,y,marco_ux,25); axis([x0,xEnd,y0,yEnd]);
            axis equal; title({['Marcoscopic Ux'];['time=',num2str(tstep)]});grid on
            subplot(2,3,6); contourf(x,y,marco_uy,25); axis([x0,xEnd,y0,yEnd]);
            axis equal; title({['Marcoscopic Uy'];['time=',num2str(tstep)]});grid on
            pause(0.5)
        end
    end
end