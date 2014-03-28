function [g,h,g_w,h_w] = boundaryf(g,h,mirco_v,weightb,nv,apha)
% 為了求解邊界問題，利用黃俊成學長之paper使得在邊界上之點達到質量、動量、能量守恆

u_w = 0;%set wall velocity is zero

%use momentum and energy conservation law get T_wall and density_wall
n = densityfunc(g,h,weightb,mirco_v,1:nv);

n_w = 2*n;
T_w = 1;

% use T_wall and density_wall and u_wall get wall distribution equation
[g_w,h_w] = f_equilibrium(u_w,mirco_v,T_w,n_w,1,1:nv);


g = apha*g_w + (1-apha)*flipdim(g,1);
h = apha*h_w + (1-apha)*flipdim(h,1);