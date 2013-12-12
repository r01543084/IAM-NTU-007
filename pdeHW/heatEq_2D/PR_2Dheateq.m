%%
%Peaceman-Rachford ADI Scheme (PR scheme) solve 2D heat equation
%by Atmosphere @ NTU 2013.12.11
%%
clear all;close all;clc;
%丁丁瞒床
dx = 0.5;
dy = 0.5;
dt = 0.01;
x0 = 0;
xEnd = 10;
y0 = 0;
yEnd = 10;
tEnd = 10;
k = 10;%heat coef.
x = x0:dx:xEnd;
y = y0:dy:yEnd;
t = 0:dt:tEnd;
u = ones(length(x),length(y));%initial condition
rx = dt*k/dx^2
ry = dt*k/dy^2
%bc
u(1,:)   = 0;%y=0放
u(end,:) = 20;%y=yEnd放

u(:,1)   = 0;%x=0放
u(:,end) = 10;%x=xEnd放
%%
%x direction Thomas Algorithm solver
ajx = ones(1,length(u(1,:))-3)*(-rx/2);
bjx = ones(1,length(u(1,:))-3)*(-rx/2);
djx = ones(1,length(u(1,:))-2)*(1+rx);

for i = 2:length(djx)
    djx(i) = djx(i)-bjx(i-1)/djx(i-1)*ajx(i-1);%get new dj
end

%y direction Thomas Algorithm solver
ajy = ones(1,length(u(1,:))-3)*(-ry/2);
bjy = ones(1,length(u(1,:))-3)*(-ry/2);
djy = ones(1,length(u(1,:))-2)*(1+ry);

for i = 2:length(djy)
    djy(i) = djy(i)-bjy(i-1)/djy(i-1)*ajy(i-1);%get new dj
end
%%
%main loop
for n = 1:length(t)
    

    %step one
    %x direction
    for i = 2:length(x)-1
        uTx = u(i,2:end-1);
        cjx = uTx;
        cjx(1) = cjx(1)+rx/2*u(i,1);
        cjx(end) = cjx(end)+rx/2*u(i,end);
        %===============solving loop==========
        
        u(i,2:end-1) = thomas(ajx,bjx,cjx,djx,uTx);
    end
    
    %y direction
    for j = 2:length(y)-1
        u(2:end-1,j) = ry/2*(u(2:end-1,j+1)-2*u(2:end-1,j)+u(2:end-1,j-1))+u(2:end-1,j);
    end
    
    %step two
    %x direction    
    for i = 2:length(x)-1
        u(i,2:end-1) = rx/2*(u(i+1,2:end-1)-2*u(i,2:end-1)+u(i-1,2:end-1))+u(i,2:end-1);
    end
    
	%x direction
    for j = 2:length(y)-1
        uTy = u(2:end-1,j);
        cjy = uTy;
        cjy(1) = cjy(1)+ry/2*u(1,j);
        cjy(end) = cjy(end)+ry/2*u(end,j);
        %===============solving loop==========
        
        u(2:end-1,j) = thomas(ajy,bjy,cjy,djy,uTy);
    end 
	[C h]=contourf(u);
    clabel(C,h)
    caxis([0 20])
    contourcbar
    dt*n
    pause()
    
end