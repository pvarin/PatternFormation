function D = nonUniformDerivative(x,Li,Lf,bc)
    N = length(x);

    % assumes even functions
    dx = [x(2)+x(1)-2*Li; x(3:end)-x(1:end-2); 2*Lf-x(end)-x(end-1)];
    
    % Neumann Boundary Conditions
    lowerDiag = 1./[-dx(2:end);0];
    upperDiag = 1./[0; dx(1:end-1)];
    
    D = spdiags([lowerDiag upperDiag], [-1 1],N,N);
    D(1,1) = -D(1,2);
    D(end,end) = -D(end,end-1);
end