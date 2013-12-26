function [z0,u0,T0] = IC(nx)
    %SOD's problem
	p       = [1    0.1  ];%pressure
	u       = [0    0    ];%velocity
	density = [1    0.125];%density
	E = 3/2*p+(0.5).*density.*u.^2; % Energy
    T = 4*E./density-2*u.^2;    % Temperature
    z = density./sqrt(pi*T);    % Fugacity
    
    %將x分兩段
    x_middle = ceil(nx/2);
    stage1 = 1:x_middle; stage2 = x_middle+1:nx;
    
    % Initial Condition for our 2D domain
    % Fugacity
    z0(stage1) = z(1); % region 1
    z0(stage2) = z(2); % region 2
    % Velovity in x
    u0(stage1) = u(1); % region 1
    u0(stage2) = u(2); % region 2
    % temperature
    T0(stage1) = T(1); % region 1
    T0(stage2) = T(2); % region 2