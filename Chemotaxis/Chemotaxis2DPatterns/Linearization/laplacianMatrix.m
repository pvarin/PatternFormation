function L = laplacianMatrix(N,dx,bc)
    E = ones(N,1);
    L = spdiags([E -2*E E],-1:1,N,N);
    
    switch upper(bc)
        case 'NEUMANN'
            L(1,1) = -1;
            L(end,end) = -1;
        case 'DIRICHLET'
            L(1,1) = -3;
            L(end,end) = -3;
        case 'PERIODIC'
            L(end,1) = 1;
            L(1,end) = 1;
    end
    
    L = L/dx^2;
end