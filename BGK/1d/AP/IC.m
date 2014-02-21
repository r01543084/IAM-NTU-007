function [u0,T0,density0,p0] = IC(nx)
    %SOD's problem
        p       = [1    0.1  ];%pressure
        u       = [0    0    ];%velocity
        density = [1    0.125];%density
    T = p./(density);    % Temperature
    
    %¨ú¤¤ÂI
    x_middle = ceil(nx/2);
    stage1 = 1:x_middle; stage2 = x_middle+1:nx;
    
    % Initial Condition for our 2D domain
    % Velovity in x
    u0(stage1) = u(1); % region 1
    u0(stage2) = u(2); % region 2
    % temperature
    T0(stage1) = T(1); % region 1
    T0(stage2) = T(2); % region 2
    % density
    density0(stage1) = density(1); % region 1
    density0(stage2) = density(2); % region 2
    % pressure
    p0(stage1) = p(1); % region 1
    p0(stage2) = p(2); % region 2
    