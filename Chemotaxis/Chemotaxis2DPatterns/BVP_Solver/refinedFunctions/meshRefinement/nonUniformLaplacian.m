function Lap = nonUniformLaplacian(x,Li,Lf,bc)
    N = length(x);

    % assumes even functions
    dx = diff(x);
    
    % Neumann Boundary Conditions
    D1 = [0; 1./dx];
    D0 = [1./dx; 0];

    lowerDiag = D0;
    centerDiag = (-D0-D1);
    upperDiag = D1;
    
    d=2./[x(2)+x(1)-2*Li; x(3:end)-x(1:end-2); 2*Lf-x(end)-x(end-1)];
    lowerDiag = lowerDiag.*[d(2:end); 0];
    centerDiag = centerDiag.*d;
    upperDiag = upperDiag.*[0;d(1:end-1)];
    
    Lap = spdiags([lowerDiag centerDiag upperDiag],-1:1,N,N);
end