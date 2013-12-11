%%
%Peaceman-Rachford ADI Scheme (PR scheme) solve 2D heat equation
%by Atmosphere @ NTU 2013.12.11
%%
clear all;close all;clc;
%丁丁瞒床
dx = 0.5;
dy = 0.5;
dt = 0.001;
x0 = 0;
xEnd = 10;
y0 = 0;
yEnd = 10;
tEnd = 10;
k = 50;%heat coef.
x = x0:dx:xEnd;
y = y0:dy:yEnd;
t = 0:dt:tEnd;
u = zeros(length(x),length(y));%initial condition
rx = dt*k/dx^2
ry = dt*k/dy^2
%bc
u(1,:)   = 0;%y=0放
u(end,:) = 20;%y=yEnd放

u(:,1)   = 0;%x=0放
u(:,end) = 10;%x=xEnd放
%%
%x direction Thomas Algorithm solver
ajx = ones(1,length(u(1,:))-3)*(rx);
bjx = ones(1,length(u(1,:))-3)*(rx);
djx = ones(1,length(u(1,:))-2)*(1+2*rx);

for i = 2:length(djx)
    djx(i) = djx(i)-bjx(i-1)/djx(i-1)*ajx(i-1);%get new dj
end

%y direction Thomas Algorithm solver
ajy = ones(1,length(u(1,:))-3)*(ry);
bjy = ones(1,length(u(1,:))-3)*(ry);
djy = ones(1,length(u(1,:))-2)*(1+2*ry);

for i = 2:length(djy)
    djy(i) = djy(i)-bjy(i-1)/djy(i-1)*ajy(i-1);%get new dj
end
%%
%main loop
for n = 1:length(t)
    
    uT = u(2:end-1,:);%ui,j
    uTp = circshift(uT,[-1 0]);%ui,j+1
    uTpp = circshift(uT,[-2 0]);%ui,j+2
    cj1 = ry*uTpp+(1+2*ry)*uTp+ry*uT;
    cj1(:,end-1:end) = [];%程ㄢぃ衡ず
    
        %===============solving loop==========
	u(2:end-1,:) = thomas2d(ajx,bjx,cj1,djx,uT);
    
%     %y direction
%     for i = 2:length(y)-1
%         umy  = circshift(u(:,i),[-1 0]);%u2
%         ummy = circshift(u(:,i),[-2 0]);%u3
%         uTy = u(2:end-1,i);
%         cjy = ry*u(:,i)+(1-2*ry)*umy+ry*ummy;
%         cjy(end-1:end) = [];%程ㄢぃ衡ず
%         cjy(1) = cjy(1)+ry*u(1,i);
%         cjy(end) = cjy(end)+ry*u(end,i);
%         %===============solving loop==========
%         
%         u(2:end-1,i) = thomas(ajy,bjy,cjy,djy,uTy);
%     end
    
	[C h]=contourf(u);
    clabel(C,h)
    caxis([0 20])
    contourcbar
    pause()
    
end