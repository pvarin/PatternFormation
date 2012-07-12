v0 = 1;
kappa = 1;
dx = .02;
mu = -3;

[u v] = full_chemIVPShoot(1,1,dx,mu);

N = length(u);
L = dx*N*2;
mass = mean(u);


guess = [u;v;0];
bvp = full_ChemBvp(kappa,N);
sol = fsolve(@(sol) bvp(sol,mass,L),guess);

plot(bvp(sol,mass,L))
title('Error in the BVP');

U = sol(1:(end-1)/2);
V = sol((end+1)/2:end-1);
sigma = sol(end);

[err A] = linearizedOperator(1,U,V,L,'Neumann',2);
[e_vect e_val] = eigs(A,5,'lr');
diag(e_val)