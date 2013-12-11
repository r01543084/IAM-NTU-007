% Solves numerically a two-dimensional diffusion equation, using forward-
% time difference, backward-time difference, Crank-Nicolson, or Peaceman-
% Rachford method.
% Francisco J. Beron-Vera, 2006/15/05
clear
disp('Solves numerically the two-dimensional diffusion equation:')
disp(' ')
disp(' u_t = u_{xx} +u_{yy} (x,y,t) \in [0,1]^2 \times [0,0.1],')
disp(' u(x,y,0) = \sin \pi x \sin \pi y (x,y) \in [0,1]^2,')
disp(' u(0,y,t) = 0 (y,t) \in [0,1] \times [0,0.1],')
disp(' u(1,y,t) = 0 (y,t) \in [0,1] \times [0,0.1],')
disp(' u(x,0,t) = 0 (x,t) \in [0,1] \times [0,0.1],')
disp(' u(x,1,t) = 0 (x,t) \in [0,1] \times [0,0.1],')
disp(' ')
disp('whose exact solution is')
disp(' ')
disp(' u(x,y,t) = \sin \pi x \sin \pi y \exp(-2\pi^2 t).')
disp(' ')
disp('Available methods:')
disp(' ')
disp('(1) Forward-time difference : U(t+Dt) = Pfd*U(t) (N.B. Dx^2 + Dy^2 > 4*Dt)')
disp('(2) Backward-time difference: Pbd*U(t+Dt) = U(t)')
disp('(3) Crank-Nicolson : Pcn1*U(t+Dt) = Pcn2*U(t)')
disp('(4) Peaceman-Rachford : Ppr1x*U(t+Dt) = Ppr2y*U(t)')
disp(' Ppr1y*U(t+Dt) = Ppr2x*U(t+Dt)')
disp(' ')
method = input('Your method choice: ');
Nx = input('Number of grid points in x direction: ') - 1;
Ny = input('Number of grid points in y direction: ') - 1;
Nt = input('Number of time steps = ');
u = @(x,y,t) sin(pi*x).*sin(pi*y)*exp(-2*pi^2*t);
Dx = 1/Nx;
x = (0:Nx)'*Dx; Dy = 1/Ny;
y = (0:Ny)'*Dy; Dt = .1/Nt;
t = (0:Nt)'*Dt; sx = Dt/Dx^2;
sy = Dt/Dy^2; nx = Nx-1; ny = Ny-1;
Sx = repmat(sx,nx,1); Sy = repmat(sy,ny,1);
Z = spalloc(nx*ny,1,nx*ny); 
E = speye(nx*ny);
Ex = speye(nx); Ey = speye(ny);
d2dx2 = spdiags([Sx -2*Sx Sx],-1:1,nx,nx);
d2dy2 = spdiags([Sy -2*Sy Sy],-1:1,ny,ny);
nabla2 = kron(d2dy2,Ex) + kron(Ey,d2dx2);
nabla2x = spdiags([Z diag(nabla2)+2*sy Z],[-nx 0 nx],nabla2);
nabla2y = spdiags([Z diag(nabla2)+2*sx Z],-1:1,nabla2);
if method == 1
    Pfd = E + nabla2;
elseif method == 2
    Pbd = E - nabla2;
elseif method == 3
    Pcn1 = E - nabla2/2;
    Pcn2 = E + nabla2/2;
else Ppr1x = E - nabla2x/2;
    Ppr2x = E + nabla2x/2; 
    Ppr1y = E - nabla2y/2;
    Ppr2y = E + nabla2y/2;
end
% Initializing U...
[x y] = ndgrid(x,y);
UdU = u(x,y,0);
subplot(221)
pcolor(x,y,UdU)
axis([0 1 0 1])
shading flat, caxis([0 1])
set(gca,'XTick',[0 .5 1],'YTick',[0 .5 1])
xlabel('x'), ylabel('y'), title('NUMERIC')
text(1.15,1.125,'t = 0','Hor','cen')
subplot(222)
pcolor(x,y,UdU)
axis([0 1 0 1])
shading flat, caxis([0 1])
h = colorbar('SouthOutside');
set(h,'Pos',[.3 .45 .4 .025])
set(get(h,'XLabel'),'Str','u')
set(gca,'XTick',[0 .5 1],'YTick',[0 .5 1])
xlabel('x'), title('EXACT')
drawnow 
U = UdU;
U([1 end],:) = [];
U(:,[1 end]) = [];
U = U(:);
% Time stepping U...
for J = 2:Nt+1
    if method == 1
        U = Pfd*U;
    elseif method == 2
        U = Pbd\U;
    elseif method == 3 
        U = Pcn1\(Pcn2*U);
    else
        U = Ppr1x\(Ppr2y*U);
        U = Ppr1y\(Ppr2x*U);
    end
    subplot(221)
    UdU(2:end-1,2:end-1) = reshape(U,[nx ny]);
    pcolor(x,y,UdU)
    axis([0 1 0 1])
    shading flat,
    caxis([0 1])
    set(gca,'XTick',[0 .5 1],'YTick',[0 .5 1])
    xlabel('x'), ylabel('y'), title('NUMERIC')
    text(1.15,1.125,['t = ' num2str(t(J))],'Hor','cen')
    subplot(222)
    pcolor(x,y,u(x,y,t(J)))
    axis([0 1 0 1])
    shading flat, caxis([0 1])
    h = colorbar('SouthOutside');
    set(h,'Pos',[.3 .45 .4 .025])
    set(get(h,'XLabel'),'Str','u')
    set(gca,'XTick',[0 .5 1],'YTick',[0 .5 1])
    xlabel('x'), title('EXACT')
    drawnow
end