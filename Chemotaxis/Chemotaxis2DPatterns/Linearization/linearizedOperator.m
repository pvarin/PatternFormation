function [err A] = linearizedOperator(kappa,V_equil,L,bc,numHalfPeriods)
    
    if numHalfPeriods == 1
        V = V_equil;
    elseif numHalfPeriods == 2
        V = [V_equil(1:end); V_equil(end:-1:1)];
    elseif numHalfPeriods == 4
        V = [V_equil(1:end); V_equil(end:-1:1); V_equil(1:end); V_equil(end:-1:1)];
    end
    % Construct the Laplacian and Derivative operators
    N = length(V);
    dx = L/(2*length(V_equil));
    
    Lap = laplacianMatrix(N,dx,bc);
    D = derivativeMatrix(N,dx,bc);
    
    % Find the derivatives of the equilibria
    V_ = V;
    U_ = exp(V_)*mean(V_)/mean(exp(V_));
    
    V_x = D*V_;
    U_x = U_.*V_x;
    
    V_xx = V_-U_;
    
    
    
    
    % Construct the linearized operator in parts
    I = eye(N);
    
    V_x = spdiags(V_x,0,N,N);
    V_xx = spdiags(V_xx,0,N,N);
    
    U_ = spdiags(U_,0,N,N);
    U_x = spdiags(U_x,0,N,N);
    
    A = [Lap-V_x*D-V_xx ,   -U_*Lap-U_x*D;
         I              ,   kappa*Lap-I];
     
    err=max(Lap*U_-D*(U_.*V_x));% kappa*Lap*V_+U_-V_
end