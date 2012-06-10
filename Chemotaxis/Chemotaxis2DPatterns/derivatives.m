function res = derivatives(t,X)
    kappa = 1;
    mu = 1;
    
    V = X(1);
    Vx = X(2);
    Vxx = -1/kappa*(V+exp(mu+V));
    
    res = [Vx; Vxx];
end