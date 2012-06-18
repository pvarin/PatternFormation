function A = linearizedOperator(kappa,V_equil,L,bc,numHalfPeriods)

    if numHalfPeriods == 2
        V_equil = [V_equil(1:end); V_equil(end:-1:1)];
    elseif numHalfPeriods == 4
        V_equil = [V_equil(1:end); V_equil(end:-1:1); V_equil(1:end); V_equil(end:-1:1)];
    end
    % Construct the Laplacian and Derivative operators
    N = length(V_equil);
    dx = L/(2*N);
    
    Lap = laplacianMatrix(N,dx,bc);
    D = derivativeMatrix(N,dx,bc);
    
    % Find the derivatives of the equilibria
    V_ = V_equil;
    V_x = D*V_;
    V_xx = Lap*V_;
    
    U_ = exp(V_)*mean(V_)/mean(exp(V_));
    U_x = U_.*V_x;
    
    % Construct the lineaized operator in parts
    I = zeros(N);
    
    V_x = spdiags(V_x,0,N,N);
    V_xx = spdiags(V_xx,0,N,N);
    
    U_ = spdiags(U_,0,N,N);
    U_x = spdiags(U_x,0,N,N);
    
    A = [L-V_x*D-V_xx   ,   -U_*L-U_x*D;
         kappa*L-I      ,   I          ];
end