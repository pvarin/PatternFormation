% model parameters
a = .5;
kappa = .1;
gamma = 1;

% similation parameters
dt = .1; T = 1100;
dx = 1;
numsteps = T/dt;
N = 100;

% initial conditions
C0 = zeros(N,1);
E0 = a*ones(N,1); E0(1,1) = E0(1,1) + .01;

% initialize the matrices
C = zeros(N,numsteps); E = zeros(N,numsteps);
C(:,1) = C0; E(:,1) = E0;

% construct calculation matrix
uDiag = [0; ones(N-1,1)];
cDiag = [-1;-2*ones(N-2,1);-1];
lDiag = [ones(N-1,1); 0];
B = [uDiag, cDiag, lDiag];
d = [1; 0; -1];
L = spdiags(B, d, N, N)*dt/(2*dx^2);

[L_C, U_C, ~] = lu(speye(N)-L);
[L_E, U_E, ~] = lu(speye(N)-kappa*L);

%Explicit Euler Method
for i=2:numsteps
    % domain interior
    f = E(:,i-1).*(1-E(:,i-1)).*(E(:,i-1)-a)+gamma*C(:,i-1);
    C(:,i) = U_C\(L_C\(C(:,i-1)-f*dt));
    E(:,i) = U_E\(L_E\(E(:,i-1)+f*dt));
end

figure
pcolor(linspace(0,dx*N,N),linspace(0,dt*numsteps,100),C(:,1:floor(end/100):end)')
shading flat
drawnow
    
animate(C,E,dx)