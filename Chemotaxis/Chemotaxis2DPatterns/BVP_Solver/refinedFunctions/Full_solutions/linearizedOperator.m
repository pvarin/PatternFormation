function [err A] = linearizedOperator(kappa,U_equil,V_equil,L,bc,numHalfPeriods)
    
    % Construct the Laplacian and Derivative operators
    N = length(V_equil);
    dx = L/(2*N);
    
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
    
    if numHalfPeriods == 1
        V_ = V_equil;
        U_ = U_equil;
        
    elseif numHalfPeriods == 2
        V_ = [V_equil(1:end); V_equil(end:-1:1)];
        U_ = [U_equil(1:end); U_equil(end:-1:1)];
        
        Dl = [Dl zeros(N);
              zeros(N) -Dl(end:-1:1,end:-1:1)];
        Dr = [Dr zeros(N);
              zeros(N) -Dr(end:-1:1,end:-1:1)];
        
    elseif numHalfPeriods == 4
        V_ = [V_equil(1:end); V_equil(end:-1:1); V_equil(1:end); V_equil(end:-1:1)];
        U_ = [U_equil(1:end); U_equil(end:-1:1); U_equil(1:end); U_equil(end:-1:1)];
        
        Dl = [Dl zeros(N,3*N);
              zeros(N) -Dl(end:-1:1,end:-1:1) zeros(N,2*N);
              zeros(N,2*N) Dl zeros(N);
              zeros(N,3*N) -Dl(end:-1:1,end:-1:1)];
        Dr = [Dr zeros(N,3*N);
              zeros(N) -Dr(end:-1:1,end:-1:1) zeros(N,2*N);
              zeros(N,2*N) Dr zeros(N);
              zeros(N,3*N) -Dr(end:-1:1,end:-1:1)];
    end
    
    Lap = laplacian(dx,length(V_),bc);
    
    
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