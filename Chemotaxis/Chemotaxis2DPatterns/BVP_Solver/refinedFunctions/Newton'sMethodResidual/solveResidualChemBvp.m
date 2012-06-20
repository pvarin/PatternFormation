%% solveResidualChemBvp.m
% uses fsolve to try to decrease the error in the equilibrium solution to
% the chemotaxis equations
function [U V sigma err] = solveResidualChemBvp(kappa,V,U,L,mass)
    
    N = length(V);
    guess = [U;V;0];
    
    % Establish the bvp to solve and call fsolve
    bvp = residualChemBvp(kappa,guess,.9);
    sol = fsolve(@(sol) bvp(sol,mass,L),guess);
    
    % Extract U V and sigma
    U = sol(1:(end-1)/2);
    V = sol((end+1)/2:end-1);
    sigma = sol(end);
    
    % Calculate the error
    err = bvp(sol,mass,L);
end