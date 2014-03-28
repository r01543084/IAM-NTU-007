function [g,h,g_w,h_w] = boundaryf(g,h,mirco_v,weightb,nv,apha)
% ���F�D����ɰ��D�A�Q�ζ��T���Ǫ���paper�ϱo�b��ɤW���I�F���q�B�ʶq�B��q�u��

u_w = 0;%set wall velocity is zero

%use momentum and energy conservation law get T_wall and density_wall
n = densityfunc(g,h,weightb,mirco_v,1:nv);

n_w = 2*n;
T_w = 1;

% use T_wall and density_wall and u_wall get wall distribution equation
[g_w,h_w] = f_equilibrium(u_w,mirco_v,T_w,n_w,1,1:nv);


g = apha*g_w + (1-apha)*flipdim(g,1);
h = apha*h_w + (1-apha)*flipdim(h,1);