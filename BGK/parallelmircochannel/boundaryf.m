function [gT,hT,gB,hB,gT_w,hT_w,gB_w,hB_w] = boundaryf(gT,hT,gB,hB,apha,lengthbx,lengthby,...
                                                 id_spacebx,id_spaceby,idvxbx,idvybx,idvxby,idvyby,...
                                                 mirco_vx,mirco_vy,weightx,weighty,weightxb,weightyb)
% 為了求解邊界問題，利用黃俊成學長之paper使得在邊界上之點達到質量、能量守恆

% boundary of x dir Velocity
ux_wbx = zeros(1,lengthby);%set wall velocity is zero
uy_wbx = zeros(1,lengthby);%set wall velocity is zero

% boundary of y dir Velocity
ux_wby = zeros(1,lengthbx);%set wall velocity is zero
uy_wby = zeros(1,lengthbx);%set wall velocity is zero

%use mass conservation law get T_wall and density_wall
nT = densityfunc(gT,hT,weightx,mirco_vx,flipdim(weightyb,1),mirco_vy,idvxbx,idvybx);
nB = densityfunc(gB,hB,weightx,mirco_vx,weightyb,mirco_vy,idvxbx,idvybx);

nT_w = 2*nT;% for wall top
nB_w = 2*nB;% for wall bottom

% wall temperature
T_wy = 1*ones(1,lengthby)*0.1;%for wall slop = NAN
T_wx = 1*ones(1,lengthbx)*0.1;%for wall slop = 0

% use T_wall and density_wall and u_wall get wall distribution equation
[gT_w,hT_w] = f_equilibrium(ux_wby,mirco_vx,uy_wby,mirco_vy,T_wx,nT_w,id_spacebx,idvxbx,idvybx);
[gB_w,hB_w] = f_equilibrium(ux_wby,mirco_vx,uy_wby,mirco_vy,T_wx,nB_w,id_spacebx,idvxbx,idvybx);

% reflaction(0)>>>>diffusitive(1), apha
% renew distribution equation
gB = apha*gB_w + (1-apha)*flipdim(gB,4);
hB = apha*hB_w + (1-apha)*flipdim(hB,4);

gT = apha*gT_w + (1-apha)*flipdim(gT,4);
hT = apha*hT_w + (1-apha)*flipdim(hT,4);

[gT_w,hT_w] = f_equilibrium(ux_wby,mirco_vx,uy_wby,mirco_vy,T_wx,nT_w,id_spacebx,idvxbx,idvybx);
[gB_w,hB_w] = f_equilibrium(ux_wby,mirco_vx,uy_wby,mirco_vy,T_wx,nB_w,id_spacebx,idvxbx,idvybx);