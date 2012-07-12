function shootTest
    kappa = 2;
    mu = -2;
    
    function dU_dx = derivatives(~,U)
        u = U(1);
        u_x = U(2);
        
        u_xx = u_x/u + u*(-u+log(u)-mu)/kappa;
        
        
        dU_dx = [u_x; u_xx];
    end

    [X,U] = ode45(@derivatives, [0 20],[exp(-1.8) 0]);
    plot(X,U)
    shg
end