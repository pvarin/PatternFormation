function err = Error(V,L,kappa)
    mu = log(mean(V)/mean(exp(V)));
    U = exp(V+mu);
    
    Lap = laplacianMatrix(length(V),L/(2*length(V)),'Neumann');
    D = derivativeMatrix(length(V),L/(2*length(V)),'Neumann');
    
    U_t = D*D*U-D*(U.*(D*V));
    U_t = Lap*U-D*U.*(D*V)-U.*(Lap*V);
    V_t = kappa*Lap*V+U-V;
    V_t = kappa*D*D*V+U-V;
    
    err = [U_t V_t];
%     err = [U_t'*U_t, V_t'*V_t];
end