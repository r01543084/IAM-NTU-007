function [z0,u0,T0] = IC(nx,ny)
    %SOD's problem
	p       = [1    0.1;
               1    0.1];%pressure
	u       = [0    0  ;    
               0    0];%velocity
	density = [1    0.125;
               1    0.125];%density
           
	E = 3/2*p+(0.5).*density.*u.^2; % Energy
    T = 4*E./density-2*u.^2;    % Temperature
    z = density./sqrt(pi*T);    % Fugacity
    
    %將x,y分兩段
    x_middle = ceil(nx/2);
    stagex1 = 1:x_middle; stagex2 = x_middle+1:nx;
	y_middle = ceil(ny/2);
    stagey1 = 1:y_middle; stagey2 = y_middle+1:ny;
    
    % Initial Condition for our 2D domain
    % Fugacity
    z0(:,stagex1,stagey1) = z(1,1); % region 1
    z0(:,stagex1,stagey2) = z(1,2); % region 2
    z0(:,stagex2,stagey1) = z(2,1); % region 3
    z0(:,stagex2,stagey2) = z(2,2); % region 4
    % Velovity in x
    u0(:,stagex1,stagey1) = u(1,1); % region 1
    u0(:,stagex1,stagey2) = u(1,2); % region 2
    u0(:,stagex2,stagey1) = u(2,1); % region 3
    u0(:,stagex2,stagey2) = u(2,2); % region 4
    % temperature
	T0(:,stagex1,stagey1) = T(1,1); % region 1
    T0(:,stagex1,stagey2) = T(1,2); % region 2
    T0(:,stagex2,stagey1) = T(2,1); % region 3
    T0(:,stagex2,stagey2) = T(2,2); % region 4
