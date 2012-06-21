function [u L N] = chemIVPShoot(u0,dx,kappa,mu)
    u(1) = u0;
    u(2) = u0;
    i=2;
    
    while u(i) >= u(i-1)
        % symmetric difference
        a = -1/(4*dx^2*u(i)^2);
        b = 1/(u(i)*dx) + u(i-1)/(2*dx^2*u(i)^2);
        c = u(i) - u(i-1)^2/(4*dx^2*u(i)) - log(u(i)) + mu;
        
        u(i+1) = (-b+sqrt(b^2-4*a*c))/(2*a);

        % left-sided derivative
%         u(i+1) = dx^2*u(i)*(log(u(i))-mu-u(i))/kappa - ...
%                   u(i)*(1-u(i-1)/u(i)) + ...
%                   2*u(i)-u(i-1);
        
        i = i+1;
        plot(u);drawnow
    end
    u = u';
    N = length(u);
    L = 2*dx*N;
end