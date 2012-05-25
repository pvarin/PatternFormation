% model parameters
b = 1.01;
kappa = 5;

% simulation parameters
dt = .001; T = 80;
dx = .1;
numsteps = T/dt;
N = 1000;

% initial conditions
U0 = ones(N,1); U0(1) = U0(1) + .1;
V0 = ones(N,1);

% initialize the matrices
U = zeros(N,numsteps); V = zeros(N,numsteps);
U(:,1) = U0; V(:,1) = V0;

%construct derivative matrix
B = [-ones(N,1) ones(N,1)];
d = [-1; 1];
D = spdiags(B, d, N, N);
D(1,1) = -1;
D(N,N) = 1;
D = D/dx;

% construct laplacian matrix
uDiag = [0; ones(N-1,1)];
cDiag = [-1;-2*ones(N-2,1);-1];
lDiag = [ones(N-1,1); 0];
B = [uDiag, cDiag, lDiag];
d = [1; 0; -1];
L = spdiags(B, d, N, N)/dx^2;

[L_U, U_U, ~] = lu(speye(N)-dt/2*L);
[L_V, U_V, ~] = lu(speye(N)-dt/2*(kappa*L-speye(N)));

U_op = speye(N)+dt/2*L;
V_op = speye(N)+dt/2*(kappa*L-speye(N));

%Explicit Euler Method
for i=2:numsteps
    U_op*U(:,i-1);
    U(:,i) = U_U\(L_U\(U_op*U(:,i-1)-dt*b*D*(D*V(:,i-1).*U(:,i-1))));
    V(:,i) = U_V\(L_V\(V_op*V(:,i-1)+dt*U(:,i-1)));
end

spacetime(U(:,1:floor(end/100):end),linspace(0,dx*N,N),linspace(0,dt*numsteps,100))   
animate(U,V,dx)

U = U(:,1:80:80000); V = V(:,1:80:80000);
x = dx:dx:N*dx; t = dt*80:dt*80:80;
save('CrankNicolsonChemotaxisData_new.mat','U','V','x','t','dx','dt','b','kappa')