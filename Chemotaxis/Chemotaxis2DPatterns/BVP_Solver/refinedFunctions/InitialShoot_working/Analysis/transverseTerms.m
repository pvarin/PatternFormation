function M = transverseTerms(u,v)
    N = length(u);
    
    I = speye(N);
    O = sparse(N,N);
    
    U = spdiags(u,0,N,N);
    V = spdiags(v,0,N,N);
    
    M = [-I, U;
         O, -I];
end