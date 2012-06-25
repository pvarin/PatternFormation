%% Create an initial guess using explicit iteration
% mu=1-e, v0=vmax=e+0.8, gives periodic solution of period approximately 6.1440 with 1024 points
% use kappa=1, dx=0.003
kappa=1; dx=0.0030008;
% kappa = 1; dx = 0.0060016;
[v N L]=initialshootfull(exp(1)+0.8,dx,mu,kappa);
u=exp(mu+v);
% u = [u; u(end:-1:1)];
% v = [v; v(end:-1:1)];
init=[u; v];
% N = N*2;
% L = 2*L;

%% Confirm this solution with chembvp
fullerr=fullchembvp(init,N,mu,L,kappa);
plot(fullerr)
max(fullerr)
mean(abs(fullerr))

%% Refine the solution using chemsolve alias fsolve
vsol=chemsolvefull(init,N,mu,L,kappa);

%% Construct the Linearization and Compute the Spectrum
dx=L/N;
Lap = laplacian(dx,N,'periodic');
D = derivative(dx,N,'symmetric','periodic');

%  rhs=Lap*v+exp(mu+v)-v;%  
u=vsol(1:N);v=vsol(N+1:2*N);

Linchem=  [Lap-spdiags(D*(D*v),0,N,N)-spdiags((D*v),0,N,N)*D  ...
      -spdiags((D*u),0,N,N)*D-spdiags(u,0,N,N)*D*D;
     speye(N,N) ...
     kappa*Lap-speye(N,N)];
[V,Di]=eigs(Linchem,10,1);
diag(Di)
plot([V(1:N,1), V(N+1:end,1)])