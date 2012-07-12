function L = laplacian(dx,N,bc)
    %% Construct the Laplacian Matrix
    E = ones(N,1);
    L = spdiags([E -2*E E],-1:1,N,N);
    L = L/dx^2;
    
    %% Employ the Boundary Conditions
    % considers boundaries that fall on the half grid (exactly dx/2 from
    % the first and last mesh points)
    switch upper(bc)
        case 'NEUMANN'
            L(1,1) = -1;
            L(end,end) = -1;
        case 'DIRICHLET'
            L(1,1) = -3;
            L(end,end) = -3;
        case 'PERIODIC'
            L(1,end) = 1;
            L(end,1) = 1;
    end
end