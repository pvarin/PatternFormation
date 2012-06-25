function [err A] = linearizedOperator_symmetric(kappa,U_equil,V_equil,L,bc,numPeriods)
    
    if numPeriods == 1
        V_ = V_equil;
        U_ = U_equil;
    elseif numPeriods == 2
        V_ = [V_equil(1:end); V_equil(end:-1:1)];
        U_ = [U_equil(1:end); U_equil(end:-1:1)];
    end
    
    % Construct the Laplacian and Derivative operators
    N = length(V_);
    dx = L/(2*length(V_equil));
    
    Lap = laplacian(dx,N,bc);
    switch upper(bc)
        case 'NEUMANN'
            Dl = derivative(dx,N,'left','');
            Dr = derivative(dx,N,'right',bc);
        case 'DIRICHLET'
            Dl = derivative(dx,N,'left',bc);
            Dr = derivative(dx,N,'right','');
        case 'PERIODIC'
            Dl = derivative(dx,N,'left',bc);
            Dr = derivative(dx,N,'right',bc);
    end
    
    % Find the derivatives of the equilibria
    V_x = Dr*V_;
    U_x = Dl*U_;
    
    V_xx = V_-U_;
    
    %Calculate the error in the U equation
    err=[Lap*U_-Dl*(U_.*V_x) kappa*Lap*V_+U_-V_];
    
    
    % Construct the linearized operator in parts
    I = speye(N);
    
    V_x = spdiags(V_x,0,N,N);
    V_xx = spdiags(V_xx,0,N,N);
    
    U_ = spdiags(U_,0,N,N);
    U_x = spdiags(U_x,0,N,N);
    
    A = [Lap-V_x*Dl-V_xx ,   -U_*Lap-U_x*Dr;
         I              ,   kappa*Lap-I];
     
%     err=Lap*diag(U_)-Dl*diag((U_.*V_x));% kappa*Lap*V_+U_-V_
end