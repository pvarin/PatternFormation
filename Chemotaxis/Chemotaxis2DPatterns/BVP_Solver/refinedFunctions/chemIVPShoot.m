function [v L N] = chemIVPShoot(v0,dx,kappa,mu)
    v(1) = v0;
    v(2) = v0;
    i=2;
    
    while v(i) >= v(i-1)
        v(i+1) = 2*v(i)-v(i-1)+dx^2/kappa*(v(i)-exp(mu+v(i)));
        i = i+1;
    end
    v = v';
    N = length(v);
    L = 2*dx*N;
end