function bvp = residualChemBvp(kappa,sol_orig,factor)

    N = (length(sol_orig)-1)/2;
    
    E = ones(N,1);
    Lap = spdiags([E -2*E E],-1:1,N,N);
    Lap(1,1) = -1; Lap(end,end) = 1;
    
    D = spdiags([-E E],[-1 1],N,N);
    D(1,1) = -1; D(end,end) = 1;
    D = D/2;
    
    U = @(sol) sol(1:(end-1)/2);
    V = @(sol) sol((end+1)/2:end-1);
    dx = @(L,N) L/(2*N);
    
    res = @(L) factor*[Lap/dx(L,N)^2*U(sol_orig) - ...
            D/dx(L,N)*(U(sol_orig).*(D*V(sol_orig)));
            
            kappa*Lap/dx(L,N)^2*V(sol_orig) + ...
            U(sol_orig) - V(sol_orig)
            
            0];
    
    bvp = @(sol,mass,L) ...
        [Lap/dx(L,N)^2*U(sol) - ...
            D/dx(L,N)*(U(sol).*(D*V(sol))) + ...
            sol(end);
            
         kappa*Lap/dx(L,N)^2*V(sol) + ...
            U(sol) - V(sol);
            
            mean(U(sol))-mass] - res(L);
        
end