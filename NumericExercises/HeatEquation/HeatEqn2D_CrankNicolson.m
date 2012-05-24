% simulation parameters
dt = .01;
dx = 1; dy = 1;
numsteps = 250;
N = 100; M = 100;
k=25;

% aribitrary initial condition
[X,Y] = meshgrid(linspace(0,dx*N,N), linspace(0,dy*M,M));
V0 = X.*Y.*sin(X/10).*sin(Y/10)/10000;

% initialize the matrix
V = Matrix2Vector(V0);

% construct the 2D Laplacian matrix for Neumann boundary
diag = zeros(N*M,1);
diag(1:N) = [-2; -3*ones(N-2,1); -2];
diag(end-N+1:end) = [-2; -3*ones(N-2,1); -2];
temp = [-3; -4*ones(N-2,1); -3];%FIXME: I think this should  be
                                %[-3; -4*ones(N-2,1); -3] but this gives
                                %the wrong result, reanalyze the laplacian
                                %matrix
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
[V_L, V_U, ~] = lu(speye(N*M) - k*dt/dx^2*L);
L_op = speye(N*M) + k*dt/dx^2*L;

figure(1)

surf(linspace(0,dx*N,N),linspace(0,dx*M,M),Vector2Matrix(V,N,M))
xlabel('x')
ylabel('y')
shading interp
axis([0 N*dx 0 M*dy -.2 .8 -.2 .8])
drawnow

%Implicit Euler Method
for i=2:numsteps
    
    V = V_U\(V_L\(L_op*V));
    
    pcolor(linspace(0,dx*N,N),linspace(0,dx*M,M),Vector2Matrix(V,N,M))
    xlabel('x')
    ylabel('y')
    shading interp
    axis([0 N*dx 0 M*dy -.2 .8 -.2 .8])
    drawnow
end