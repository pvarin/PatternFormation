function rhs=fullchembvp(w,N,mu,L,kappa)

    %defines boundary value problem as function in N points using finite differences
    % u_xx - (uv_x)_x=0
    % kappa v_xx+ u-v=0

    % note that we are thinking of a periodic problem on length L with N points

    dx=L/N;
    Lap = laplacian(dx,N,'periodic');
    D = derivative(dx,N,'symmetric','periodic');

    %  rhs=Lap*v+exp(mu+v)-v;%  
    u=w(1:N);v=w(N+1:2*N);

    rhs=[(Lap*u-D*(u.*(D*v))); (kappa*Lap*v+u-v)];

    %  rhs=(Lap*u-D*(u.*(D*v)))' ;
    %  rhs=kappa*Lap*v+u-v;
end


