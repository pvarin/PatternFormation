function wsol=chemsolvefull(initial,N,mu,L,kappa)
    % finds the solutions to the boundary value problem defined in chembvp using fsolve
    % given an initial guess initial, number of points N, parameter and length L
    % 

    options = optimset('TolX',1e-4);
    wsol=fsolve(@(w) fullchembvp(w,N,mu,L,kappa),initial,options);
end