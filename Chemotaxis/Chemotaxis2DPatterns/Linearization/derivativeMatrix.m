function D = derivativeMatrix(N,dx,bc)
    E = ones(N,1);
    D = spdiags([-E E],[-1 1],N,N);
    
    switch upper(bc)
        case 'NEUMANN'
            D(1,1) = -1;
            D(end,end) = 1;
        case 'DIRICHLET'
            D(1,1) = 1;
            D(end,end) = -1;
        case 'PERIODIC'
            D(1,end) = -1;
            D(end,1) = 1;
    end
    
    D = D/(2*dx);
end