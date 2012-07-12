function [period mass] = IVP_onePeriod(kappa, b, mu, V0)
    
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

    odeOpts=odeset('Events',@events,'RelTol',1e-6);
    
    [X1 Y1] = ode45(@derivatives,[0,inf],[V0,0],opts);
    if X1(end)==inf
        period = Nan;
        mass = Nan;
        return
    end
    [X2 Y2] = ode45(@derivatives,[X1(end),inf],[Y1(end,1),Y1(end,2)],odeOpts);
    if X2(end)==inf
        period = Nan;
        mass = Nan;
        return
    end
    
    Y = [Y1; Y2];
    X = [X1; X2];
    V = Y(:,1);
    Vx = Y(:,2);
    U = exp(b*V+mu);
    H = Vx.^2/2 - (V.^2/2 - exp(mu+b*V)/b)/kappa;
    
    
%     plot(X,[V,U]);
%     ylabel('Solution');
%     xlabel('X');
%     legend('V','U')
%     figure
%     plot(X,H)
%     title('Energy of V')
%     ylabel('Energy')
%     xlabel('X')
    
    period = X(end);
    mass = trapz(X,V)/period;
end