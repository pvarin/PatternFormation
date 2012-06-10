function T = chemoPeriod(mu,kappa,b,V0)
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
    
    [X1 Y1] = ode45(@derivatives,[0,inf],[V0,0],opts);
    [X2 ~] = ode45(@derivatives,[X1(end),inf],[Y1(end,1),0],opts);
    T = X2(end);
end