function [err A] = linearizedOperator_symmetric_periodic(kappa,U_equil,V_equil,L,numPeriods)
    
    V_ = [];
    U_ = [];
    for i=1:floor(numPeriods/2)
        V_ = [V_; V_equil; V_equil(end:-1:1)];
        U_ = [U_: U_equil; U_equil(end:-1:1)];
    end
    
    if mod(numPeriods,2)% if numPeriods is odd
        V_ = [V_; V_equil];
        U_ = [U_; U_equil];
    end
    
    % Construct the Laplacian and Derivative operators
    N = length(V_);
    dx = L/(2*length(V_equil));
    
    D = derivative(dx,N,'Periodic');
    Lap = D^2;
    
    % Find the derivatives of the equilibria
    V_x = D*V_;
    U_x = D*U_;
    
    V_xx = Lap*V_xx;
    
    %Calculate the error in the U equation
    err=[Lap*U_-D*(U_.*V_x) kappa*Lap*V_+U_-V_];
    
    
    % Construct the linearized operator in parts
    I = speye(N);
    
    V_x = spdiags(V_x,0,N,N);
    V_xx = spdiags(V_xx,0,N,N);
    
    U_ = spdiags(U_,0,N,N);
    U_x = spdiags(U_x,0,N,N);
    
    A = [Lap-V_x*D-V_xx ,   -U_*Lap-U_x*D;
         I              ,   kappa*Lap-I];
     
%     err=Lap*diag(U_)-Dl*diag((U_.*V_x));% kappa*Lap*V_+U_-V_
end