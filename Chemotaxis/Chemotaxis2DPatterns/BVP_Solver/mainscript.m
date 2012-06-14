% find an initial guess using explicit iteration
% mu=1-e, v0=vmax=e+0.8, gives peruiodic solution of period approximately 6.1440 with 1024 points
% use kappa=1, dx=0.003

[v N L]=initialshoot(exp(1)+0.8,0.003,1-exp(1),1);

% confirm this solution using bvp solver in chemsolve alias fsolve
% same parameter values

vsol=chemsolve(v,N,1-exp(1),L,1);