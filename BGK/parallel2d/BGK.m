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
time_ic = [0.2 0.2 0.3 0.25 0.23 0.3 0.25 0.25 0.3 0.15 ... 
           0.3 0.25 0.3 0.1  0.2  0.2 0.3  0.2  0.3 1      ];
for ic_case = 1
%to avoid all_parameter.mat record g and h. Because they are huge.
clear g h g_star h_star g0 h0

tic
%normal coef
    name      ='CBGK2d';  % Simulation Name
    x0     = 0                  ;% X��l��m
    xEnd   = 1                  ;% X������m
    y0     = 0                  ;% Y��l��m
    yEnd   = 1                  ;% Y������m
    dx     = 1/50               ;% �Cdx���@��
    dy     = 1/50               ;% �Cdy���@��
    tEnd   = time_ic(ic_case)   ;% �q0�}�l�p��tEnd��
    r_time = 10^-8              ;% Relaxation time
    cfl    = 1                  ;% CFL nuber
    write_ans = 1               ;%(0)no,(1)yes
    draw = 2                    ;%(0)no,(1)in time, (2)last time
    using_time = 0              ;%how many time we use in this case
%plot coef.
    pic_num = 40                ;%how many picture u take
    lagtime = 0.1               ;%plot lag time
    cline = 40                  ;%contour line number
%% �����ɶ��B�t�תŶ��B��m�Ŷ�
% ��m�Ŷ�����
x = x0:dx:xEnd;
nx = length(x);

y = y0:dy:yEnd;
ny = length(y);

%�t�תŶ�����(Velocity-Space)

    %x velocity
    nvx = 20;%�]�� Gauss-Hermite ��nv���I�A���F�n���t��domain
            %order �� 2*nv-1
    [mirco_vx,weightx] = GaussHermite(nvx);%for integrating range: -inf to inf
    weightx = weightx.*exp(mirco_vx.^2);%real weight if not, chack out website
    % http://www.efunda.com/math/num_integration/findgausshermite.cfm

    %y velocity
    nvy = 20;
    [mirco_vy,weighty] = GaussHermite(nvy);
    weighty = weighty.*exp(mirco_vy.^2);

%�ɶ�����
    dt = min(dx*cfl/max(abs(mirco_vx)),dy*cfl/max(abs(mirco_vy)));
    time = [0:dt:tEnd-dt tEnd];%�j���̫�@�B��tEnd
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
        all_parameter = [ID,'/all_parameter.mat'];
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
disp(sum(sum(abs(density-density1)))*dx*dy)
disp(sum(sum(abs(marco_ux-marco_ux1)))*dx*dy)
disp(sum(sum(abs(marco_uy-marco_uy1)))*dx*dy)
disp(sum(sum(abs(T-T1)))*dx*dy)

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
        g_star(i,:,:,:) = g(i,:,:,:) + dfdx(g(i,:,:,:),mirco_vx,idvx(i,:,:,:)...
                                            ,nvx,dx,dt,[0,1])*dt;
        h_star(i,:,:,:) = h(i,:,:,:) + dfdx(h(i,:,:,:),mirco_vx,idvx(i,:,:,:)...
                                            ,nvx,dx,dt,[0,1])*dt;
                        
    end
    
    %y dir
    disp('AP. RK y dir.')
    parfor j = 1:nx
        g_star(:,j,:,:) = g_star(:,j,:,:) + dfdx(g(:,j,:,:),mirco_vy,idvy(:,j,:,:)...
                                            ,nvy,dy,dt,[1,0])*dt;
        h_star(:,j,:,:) = h_star(:,j,:,:) + dfdx(h(:,j,:,:),mirco_vy,idvy(:,j,:,:)...
                                            ,nvy,dy,dt,[1,0])*dt;
    end    
    disp('AP. RK Successful')
    
    %�Q�ηs�o�쪺f(�������)�D�o�i�o�ƱK��(n)�B�q�qor�ʶq�K��(j_x or nu)�B��q�K��
    parfor i=1:nx
        [density(:,i),marco_ux(:,i),marco_uy(:,i),T(:,i),e(:,i),p(:,i)]...
            = densityfunc(g_star(:,i,:,:),h_star(:,i,:,:),weightx,mirco_vx...
            ,weighty,mirco_vy,idvx(:,i,:,:),idvy(1:end-1,i,:,:));
    end
    
%     %���Ϯg�]Total Reflection boundary�^
%         marco_ux(21:30,21:30) = 0;
%         marco_uy(21:30,21:30) = 0;
%         density(21:30,21:30) = 0.1;
%         T(21:30,21:30) = 0.2055;
%     for i = 21:30 
%         if     marco_ux(i,20) > 10^-9	%block left side
%             marco_ux(i,20) = -marco_ux(i,20);
%         elseif marco_ux(i,31) < 10^-9   %block right side
%             marco_ux(i,31) = -marco_ux(i,31);
%         elseif marco_uy(31,i) > 10^-9   %block bottom side
%             marco_uy(31,i) = -marco_uy(31,i);
%         elseif marco_uy(i,20) < 10^-9   %block bottom side
%             marco_uy(i,20) = -marco_uy(i,20);
%         end
%     end
%     
%     if     marco_ux(20,20) > 10^-9 && marco_uy(20,20) < 10^-9%���W
%         marco_ux(20,20) = -marco_ux(20,20);
%         marco_uy(20,20) = -marco_uy(20,20);
%     elseif marco_ux(31,20) > 10^-9 && marco_uy(31,20) > 10^-9%���U
%         marco_ux(31,20) = -marco_ux(31,20);
%         marco_uy(31,20) = -marco_uy(31,20);
%     elseif marco_ux(20,31) < 10^-9 && marco_uy(20,31) < 10^-9%�k�W
%         marco_ux(20,31) = -marco_ux(20,31);
%         marco_uy(20,31) = -marco_uy(20,31);
%     elseif marco_ux(31,31) < 10^-9 && marco_uy(31,31) > 10^-9%�k�U
%         marco_ux(31,31) = -marco_ux(31,31);
%         marco_uy(31,31) = -marco_uy(31,31);
%     end
    
    %�b�C�Ӯɶ��B���A�Q�ΦU�L�[�q�A��s���źA�������
    parfor i = 1:nx
        [g_eq(:,i,:,:),h_eq(:,i,:,:)] = f_equilibrium(marco_ux...
            ,mirco_vx,marco_uy,mirco_vy,T,density,id_space(:,i,:,:)...
            ,idvx(:,i,:,:),idvy(1:end-1,i,:,:));%���źA��{��
    end
    
    %main loop use RK 4th
    
    %RK final step, max it
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
        savefunc(T,density,p,e,marco_ux,marco_uy,tstep,counter,ID)
    end
    
    if draw == 1
        %plot
        contourf_func(x,y,tstep,T,density,p,e,marco_ux,marco_uy,0.01,55);
    end
end
    
%saving this case using time
using_time = toc/3600;
tt = [ID,'/using_time','.mat'];
save(tt,'using_time','counter');

%% plot part
    % if you don't have any initial condition, then double click
    % all_parameter. And if no drawing, recheck 'draw' and 'counter' number.
    % And change ur '/Users/Atmosphere/IAM NTU 007/BGK/Copy_of_2d/' to your
    % flie add..
    if draw == 2
        for i =1:counter
            %tell matlab go where to find data
            TT=[ID,'/T/T',num2str(i),'.mat'];
            DD=[ID,'/density/density',num2str(i),'.mat'];
            pp=[ID,'/p/p',num2str(i),'.mat'];
            ee=[ID,'/e/e',num2str(i),'.mat'];
            UUx=[ID,'/Ux/Ux',num2str(i),'.mat'];
            UUy=[ID,'/Uy/Uy',num2str(i),'.mat'];
            
            %load data
            load(TT);load(DD);load(pp);load(ee);load(UUx);load(UUy);
            
            %plot
            contourf_func(x,y,tstep,T,density,p,e,marco_ux,...
                          marco_uy,lagtime,cline);
                      
            %streamline
%             [xx,yy] = meshgrid(0:1/50:1,0:1/50:1);
%             [sx,sy] = meshgrid(0:1/10:1,0:1/10:1);
%             streamline(stream2(xx,yy,marco_ux,marco_uy,sx,sy));
        end
    end
end