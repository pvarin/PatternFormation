function D = derivative(dx,N,style,bc)
    %% Construct the Derivative Matrix
    E = ones(N,1);
    
    upperDiag = zeros(N,1);
    centerDiag = zeros(N,1);
    lowerDiag = zeros(N,1);
    
    switch upper(style)
        case 'LEFT'
            centerDiag = E;
            lowerDiag = -E;
        case 'RIGHT'
            upperDiag = E;
            centerDiag = -E;
        case 'SYMMETRIC'
            upperDiag = E/2;
            lowerDiag = -E/2;
        otherwise
            warning('I do not recognize this type of derivative, please choose: "left", "right", or "symmetric"');
    end
    D = spdiags([lowerDiag centerDiag upperDiag],-1:1,N,N);
    D = D/dx;
    
    %% Employ the Boundary Conditions
    % considers boundaries that fall on the half grid (exactly dx/2 from
    % the first and last mesh points)
    switch upper(bc)
        case 'NEUMANN'
            D(1,1) = D(1,1) + lowerDiag(1);
            D(end,end) = D(end,end) + upperDiag(end);
        case 'DIRICHLET'
            D(1,1) = D(1,1) - lowerDiag(1);
            D(end,end) = D(end,end) - upperDiag(end);
        case 'PERIODIC'
            D(1,end) = D(1,end) + lowerDiag(1);
            D(end,1) = D(end,1) + upperDiag(end);
    end
end