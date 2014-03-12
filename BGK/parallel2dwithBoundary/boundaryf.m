function [gL,hL,gR,hR,gT,hT,gB,hB,g_w,h_w] = boundaryf(gL,hL,gR,hR,gT,hT,gB,hB,apha,lengthbx,lengthby,...
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
nL = densityfunc(gL,hL,weightxb,mirco_vx,weighty,mirco_vy,idvxby,idvyby);
nR = densityfunc(gR,hR,flipdim(weightxb,1),mirco_vx,weighty,mirco_vy,idvxby,idvyby);
nT = densityfunc(gT,hT,weightx,mirco_vx,flipdim(weightyb,1),mirco_vy,idvxbx,idvybx);
nB = densityfunc(gB,hB,weightx,mirco_vx,weightyb,mirco_vy,idvxbx,idvybx);

nL_w = 2*nL';% for wall left
nR_w = 2*nR';% for wall right
nT_w = 2*nT;% for wall top
nB_w = 2*nB;% for wall bottom

n_w  = [((nL_w(1)^2+nR_w(1)^2+nT_w(1))/3)^0.5 ...
        ((nL_w(2:end-1).^2+nR_w(2:end-1).^2)./2).^0.5...
        ((nL_w(end)^2+nR_w(end)^2+nB_w(1))/3)^0.5];

T_wy = 1*ones(1,lengthby);%for wall slop = NAN
T_wx = 1*ones(1,lengthbx);%for wall slop = 0

% use T_wall and density_wall and u_wall get wall distribution equation
[gL_w,hL_w] = f_equilibrium(ux_wbx,mirco_vx,uy_wbx,mirco_vy,T_wy,nL_w,id_spaceby,idvxby,idvyby);
[gR_w,hR_w] = f_equilibrium(ux_wbx,mirco_vx,uy_wbx,mirco_vy,T_wy,nR_w,id_spaceby,idvxby,idvyby);
[gT_w,hT_w] = f_equilibrium(ux_wby,mirco_vx,uy_wby,mirco_vy,T_wx,nT_w,id_spacebx,idvxbx,idvybx);
[gB_w,hB_w] = f_equilibrium(ux_wby,mirco_vx,uy_wby,mirco_vy,T_wx,nB_w,id_spacebx,idvxbx,idvybx);

% reflaction(0)>>>>diffusitive(1), apha
% renew distribution equation
gL = apha*gL_w + (1-apha)*flipdim(gL,3);
hL = apha*hL_w + (1-apha)*flipdim(hL,3);

gR = apha*gR_w + (1-apha)*flipdim(gR,3);
hR = apha*hR_w + (1-apha)*flipdim(hR,3);

gB = apha*gB_w + (1-apha)*flipdim(gB,4);
hB = apha*hB_w + (1-apha)*flipdim(hB,4);

gT = apha*gT_w + (1-apha)*flipdim(gT,4);
hT = apha*hT_w + (1-apha)*flipdim(hT,4);

[g_w,h_w] = f_equilibrium(ux_wbx,mirco_vx,uy_wbx,mirco_vy,T_wy,n_w,id_spaceby,idvxby,idvyby);
