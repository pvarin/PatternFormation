% simulation parameters
dt = .01;
dx = 1; dy = 1;
numsteps = 1000;
N = 100; M = 100;
k=10;

% aribitrary initial condition
[X,Y] = meshgrid(linspace(0,dx*N,N), linspace(0,dy*M,M));
U0 = sin(X/10).*sin(Y/10);

% initialize the matrices
U = zeros(N,M,numsteps);

U(:,:,1) = U0;

figure(1)

%Explicit Euler Method
for i=2:numsteps
    % domain interior
    n = 2:N-1;
    m = 2:M-1;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)+U(n,m+1,i-1)-4*U(n,m,i-1)+U(n-1,m,i-1)+U(n,m-1,i-1))*dt/dx^2;
    
    % domain boundaries (Neumann Conditions)
    % Edges
    n = 1;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)+U(n,m+1,i-1)-3*U(n,m,i-1)+U(n,m-1,i-1))*dt/dx^2;
    
    n = N;
    U(n,m,i) = U(n,m,i-1)+k*(U(n,m+1,i-1)-3*U(n,m,i-1)+U(n-1,m,i-1)+U(n,m-1,i-1))*dt/dx^2;
    
    n = 2:N-1;
    m = 1;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)+U(n,m+1,i-1)-3*U(n,m,i-1)+U(n-1,m,i-1))*dt/dx^2;

    m = M;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)-3*U(n,m,i-1)+U(n-1,m,i-1)+U(n,m-1,i-1))*dt/dx^2;
    
    % Corners
    n = N;
    U(n,m,i) = U(n,m,i-1)+k*(-2*U(n,m,i-1)+U(n-1,m,i-1)+U(n,m-1,i-1))*dt/dx^2;

    n = 1;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)-2*U(n,m,i-1)+U(n,m-1,i-1))*dt/dx^2;

    m = 1;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)+U(n,m+1,i-1)-2*U(n,m,i-1))*dt/dx^2;
    
    m = M;
    U(n,m,i) = U(n,m,i-1)+k*(U(n+1,m,i-1)-2*U(n,m,i-1)+U(n,m-1,i-1))*dt/dx^2;

    surf(linspace(0,dx*N,N),linspace(0,dx*M,M),U(:,:,i))
    xlabel('x')
    ylabel('y')
    shading interp
    axis([0 N*dx 0 M*dy -1 1 -1 1])
    drawnow
end