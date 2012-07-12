function [err A] = linearizedOperator_periodic(kappa,U_equil,V_equil,L,numPeriods)
    
    % Construct the Laplacian and Derivative operators
    N = length(V_equil);
    dx = L/(2*N);
    
    Dl = compoundDerivative(dx,N,'left',numPeriods);
    Dr = compoundDerivative(dx,N,'right',numPeriods);
    
    if numPeriods == 1
        V_ = [V_equil(1:end); V_equil(end:-1:1)];
        U_ = [U_equil(1:end); U_equil(end:-1:1)];
        
    elseif numPeriods == 2
        V_ = [V_equil(1:end); V_equil(end:-1:1); V_equil(1:end); V_equil(end:-1:1)];
        U_ = [U_equil(1:end); U_equil(end:-1:1); U_equil(1:end); U_equil(end:-1:1)];
    end
    
    Lap = laplacian(dx,length(V_),'periodic');
    
    
    % Find the derivatives of the equilibria    
    V_x = Dr*V_;
    U_x = Dl*U_;
    
    V_xx = V_-U_;
    
    %Calculate the error in the U equation
    err=[Dl*Dr*U_-Dl*(U_.*V_x) kappa*Lap*V_+U_-V_];
    
    
    % Construct the linearized operator in parts
    I = speye(N*2);
    
    V_x = spdiags(V_x,0,N*2,N*2);
    V_xx = spdiags(V_xx,0,N*2,N*2);
    
    U_ = spdiags(U_,0,N*2,N*2);
    U_x = spdiags(U_x,0,N*2,N*2);
    
%     A = [Lap-V_x*Dl-V_xx ,   -U_*Lap-U_x*Dr;
%          I              ,   kappa*Lap-I];
     
    A = [Dl*Dr-V_x*Dl-V_xx ,   -U_*Lap-U_x*Dr;
         I              ,   kappa*Lap-I];
end