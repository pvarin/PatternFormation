function [U V sigma] = solveFullChemBvp(kappa,V,L)
    N = length(V);
    mass = mean(V);
    U = exp(V)*mass/mean(exp(V));
    guess = [U;V;0];
    
    bvp = fullChemBvp(kappa,N);
    sol = fsolve(@(sol) bvp(sol,mass,L),guess);
    
    U = sol(1:(end-1)/2);
    V = sol((end+1)/2:end-1);
    sigma = sol(end);
end