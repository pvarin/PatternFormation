function IVP_onePeriod()
    mu=-1.1;
    b=1;
    kappa=1;
    
    function res = derivatives(~,X)
        V = X(1);
        Vx = X(2);
        Vxx = 1/kappa*(V-exp(mu+b*V));

        res = [Vx; Vxx];
    end

    function [val,isterminal,direction] = events(~,X)
        val = X(2);
        isterminal = 1;
        direction = 0;
    end

    opts=odeset('Events',@events);
    
    [X1 Y1] = ode45(@derivatives,[0,inf],[1.6,0],opts);
    [X2 Y2] = ode45(@derivatives,[X1(end),inf],[Y1(end,1),0],opts);
    
    Y = [Y1; Y2];
    X = [X1; X2];
    V = Y(:,1);
    U = exp(b*V+mu);
    plot(X,[V,U])
end