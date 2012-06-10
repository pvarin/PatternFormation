function equilibriumSolutions
    mu=1;
    b=1;
    kappa=1;

    
    function res = derivatives(t,X)
        V = X(1);
        Vx = X(2);
        Vxx = -1/kappa*(V+exp(mu+b*V));

        res = [Vx; Vxx];
    end



    [X Y] = ode45(@derivatives,[0,10],[1,0]);
    V = Y(:,1);
    U = exp(b*V+mu);
    plot(X,[V,U])
    
end