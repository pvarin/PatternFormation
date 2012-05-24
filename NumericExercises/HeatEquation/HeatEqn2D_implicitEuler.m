% simulation parameters
dt = .01;
dx = 1; dy = 1;
numsteps = 250;
N = 100; M = 100;
k=25;

% aribitrary initial condition
[X,Y] = meshgrid(linspace(0,dx*N,N), linspace(0,dy*M,M));
U0 = X.*Y.*sin(X/10).*sin(Y/10)/10000;

% initialize the matrix
U = Matrix2Vector(U0);

% construct the 2D Laplacian matrix

figure(1)

%Implicit Euler Method
for i=2:numsteps
    
    
    
    surf(linspace(0,dx*N,N),linspace(0,dx*M,M),Vector2Matrix(U,N,M))
    xlabel('x')
    ylabel('y')
    shading interp
    axis([0 N*dx 0 M*dy -.2 .8 -.2 .8])
    drawnow
end