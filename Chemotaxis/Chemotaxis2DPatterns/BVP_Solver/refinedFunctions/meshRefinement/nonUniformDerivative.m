function D = nonUniformDerivative(x,Li,Lf,bc)
    N = length(x);

%     % assumes even functions
%     dx = [x(2)+x(1)-2*Li; x(3:end)-x(1:end-2); 2*Lf-x(end)-x(end-1)];
%     
%     % Neumann Boundary Conditions
%     lowerDiag = 1./[-dx(2:end);0];
%     upperDiag = 1./[0; dx(1:end-1)];
%     
    
    
    h = [2*(x(1)-Li); x(2:end)-x(1:end-1); 2*(Lf-x(end))];
    
    % first order terms
    upperDiag = 1./(h(2:end)+h(1:end-1));
    centerDiag = zeros(length(x),1);
    lowerDiag = -1./(h(2:end)+h(1:end-1));
    
    % second order terms
    upperDiag = upperDiag + 2*(h(2:end)-h(1:end-1))./ ...
                            (h(2:end)+h(1:end-1))./ ...
                            h(2:end);
                        
    centerDiag = centerDiag + -2*(h(2:end)-h(1:end-1))./ ...
                              (h(2:end)+h(1:end-1)).* ...
                              (1./h(2:end)+1./(h(1:end-1)));
    
    lowerDiag = lowerDiag + 2*(h(2:end)-h(1:end-1))./ ...
                            (h(2:end)+h(1:end-1))./ ...
                            h(1:end-1);
    
    % Restructure Diagonals and Save Boundary Data
    initialBoundaryValue = lowerDiag(1);
    finalBoundaryValue = upperDiag(end);
    
    upperDiag = [0; upperDiag(1:end-1)];
    lowerDiag = [lowerDiag(2:end); 0];
    
    % Construct the Matrix
    D = spdiags([lowerDiag centerDiag upperDiag], -1:1,N,N);
    
    switch upper(bc)
        case 'NEUMANN'
            % Homogeneous Neumann Boundary Conditions
            D(1,1) = D(1,1) + initialBoundaryValue;
            D(end,end) = D(end,end) + finalBoundaryValue;
        case 'DIRICHLET'
            % Homogeneous Dirichlet Boundary Conditions
            D(1,1) = D(1,1) - initialBoundaryValue;
            D(end,end) = D(end,end) - finalBoundaryValue;
        case 'PERIODIC'
            % Periodic Boundary Conditions
            D(1,end) = D(1,end) + initialBoundaryValue;
            D(end,1) = D(end,1) + finalBoundaryValue;
    end
end