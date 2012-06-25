% Model Parameters
kappa = 1;
mu = 1-exp(1);
v0 = exp(1)+.8;

% Simulation Parameters
dx = 0.0030008;

% Create the initial guess
[v N L] = chemIVPShoot(v0, dx, mu, kappa);
u = exp(mu+v);
sol = [u;v];

% Confirm that this solution is close
bvp = chemBVP(N);
err = bvp(sol,L,kappa);
plot(err); drawnow
fprintf('The mean square error is\t%f\n',err'*err/N)
fprintf('The maximum error is\t\t%f\n',max(abs(err)))

% Refine the solution using fsolve
sol = chemSolve(sol,L,kappa,bvp);

%% Reassess the error
err = bvp(sol,L,kappa);
plot(err); drawnow
fprintf('The mean square error is\t%f\n',err'*err/N)
fprintf('The maximum error is\t\t%f\n',max(abs(err)))

% Linearize and compute the spectrum of the linearization
u = sol(1:N);
v = sol(N+1:end);
Lin = chemLinear(u,v,L,kappa);
[eigVect eigVal] = eigs(Lin,10,1);
diag(eigVal)