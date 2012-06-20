function bvp = fullChemBvp(kappa,N)
    E = ones(N,1);
    Lap = spdiags([E -2*E E],-1:1,N,N);
    Lap(1,1) = -1; Lap(end,end) = 1;
    
    D = spdiags([-E E],[-1 1],N,N);
    D(1,1) = -1; D(end,end) = 1;
    D = D/2;
    
    U = @(sol) sol(1:(end-1)/2);
    V = @(sol) sol((end+1)/2:end-1);
    dx = @(L,N) L/(2*N);
    
    bvp = @(sol,mass,L) ...
        [Lap/dx(L,N)^2*U(sol) - ...
            D/dx(L,N)*(U(sol).*(D*V(sol))) + ...
            sol(end);
            
         kappa*Lap/dx(L,N)^2*V(sol) + ...
            U(sol) - V(sol);
            
         mean(U(sol)) - mass];
end