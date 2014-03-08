function [g,h,g_w,h_w] = boundaryf(g,h,weightx,mirco_vx,weighty,mirco_vy,...
                                   id_space,idvx,idvy,apha,lengthby)
% 為了求解邊界問題，利用黃俊成學長之paper使得在邊界上之點達到質量、能量守恆

u_wx = zeros(1,lengthby);%set wall velocity is zero
u_wy = zeros(1,lengthby);%set wall velocity is zero

%use momentum and energy conservation law get T_wall and density_wall
n= densityfunc(g,h,weightx,mirco_vx,weighty,...
                                        mirco_vy,idvx,idvy);

n_w = 2*n;
T_w = 1*ones(1,lengthby);

% use T_wall and density_wall and u_wall get wall distribution equation
[g_w,h_w] = f_equilibrium(u_wx,mirco_vx,u_wy,mirco_vy,T_w...
                         ,n_w,id_space,idvx,idvy);

% reflaction(0)>>>>diffusitive(1), apha
g = apha*g_w + (1-apha)*flipdim(g,3);
h = apha*h_w + (1-apha)*flipdim(h,3);