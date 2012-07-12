% model parameters
kappa = .01;
a=.5;
gamma = .1;

% simulation parameters
dt = .1;
dx = .1; dy = .1;
numsteps = 5000;
N = 1000; M = 1000;

% establish initial conditions
C0 = zeros(N,M);
E0 = a*ones(N,M);
E0(floor(end/2),floor(end/2)) = E0(floor(end/2),floor(end/2))+.1;

% initialize the matrix
C = Matrix2Vector(C0);
E = Matrix2Vector(E0);

% construct the 2D Laplacian matrix for Neumann boundary
diag = zeros(N*M,1);
diag(1:N) = [-2; -3*ones(N-2,1); -2];
diag(end-N+1:end) = [-2; -3*ones(N-2,1); -2];
temp = [-3; -4*ones(N-2,1); -3];
for i=2:M-1
    diag(((i-1)*N+1):i*N) = temp;
end
B = [ones(N*M,1) ones(N*M,1) diag ones(N*M,1) ones(N*M,1)];
d = [-N -1 0 1 N];
L=spdiags(B,d,N*M,N*M);
for i=1:M-1
   L(i*N+1,i*N)=0;
   L(i*N,i*N+1)=0;
end
% LU decomposition
[C_L, C_U, ~] = lu(speye(N*M) - dt/dx^2*L);
[E_L, E_U, ~] = lu(speye(N*M) - kappa*dt/dx^2*L);

figure(1)
pcolor(linspace(0,dx*N,N),linspace(0,dx*M,M),Vector2Matrix(E,N,M))
xlabel('x')
ylabel('y')
shading interp
axis([0 N*dx 0 M*dy])
drawnow

%Implicit Euler Method
for i=2:numsteps
    f = E.*(1-E).*(E-a) + gamma*C;
    
    C = C_U\(C_L\(C-dt*f));
    E = E_U\(E_L\(E+dt*f));
    
    figure(1)
    pcolor(linspace(0,dx*N,N),linspace(0,dx*M,M),Vector2Matrix(E,N,M))
    title('E')
    xlabel('x')
    ylabel('y')
    shading interp
    colormap gray
    axis([0 N*dx 0 M*dy 0 .0001 0 1])
    drawnow
end