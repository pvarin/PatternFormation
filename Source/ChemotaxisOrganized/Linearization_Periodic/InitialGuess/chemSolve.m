%% chemSolve.m
%
% Solves the equilibrium chemotaxis boundary value problem defined in
% chemBVP.m using fsolve
%
% Requires an initial guess for the solution as well as the size of the
% period, L, and the parameter kappa

function sol=chemSolve(initial,L,kappa,bvp)
    %% Set options and call fsolve
    options = optimset('TolX',1e-4);
    sol=fsolve(@(sol_) bvp(sol_,L,kappa),initial,options);
end