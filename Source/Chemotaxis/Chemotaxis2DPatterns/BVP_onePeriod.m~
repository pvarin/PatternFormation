function BVP_onePeriod
    kappa = 1;
    mu = -4;
    b = 1;
    T = 40;
    
    function res = odefun(~,X)
        V = X(1);
        Vx = X(2);
        Vxx = 1/kappa*(V-exp(mu+b*V));

        res = [Vx; Vxx];
    end

    function res = bcfun(Xa,Xb)
        res = [Xa(2)-Xb(2)];
    end

    x = linspace(0,T);
    y = @(x_) [cos(pi*x_/T)+1.1,-pi/T*sin(pi*x_/T)];
    solinit = bvpinit(x,y);
    sol = bvp5c(@odefun, @bcfun, solinit);
    
    plot(sol.x,sol.y)
end