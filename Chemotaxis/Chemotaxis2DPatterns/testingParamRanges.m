b = 1;
kappa = 1;
chi = 1;
mu = -1.1;
%find the upper bound on mu
muMax = -1-log(b);
%find the bounds on V0 (i.e. the equilibria of the V ode)
Vequil = @(V) exp(mu+b*V)-V;
Vmin = fzero(Vequil,[0,1/b]);
Vmax = fzero(Vequil,[1/b,inf]);