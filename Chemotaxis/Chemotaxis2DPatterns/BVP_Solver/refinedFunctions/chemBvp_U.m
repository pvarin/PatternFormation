function bvp=chemBvp_U(kappa,N)
    % Defines boundary value problem as function in N points using finite differences
    % kappa*(u_xx/u - u_x/u^2) + u - ln(u) + mu=0

    % We are solving for even periodic solutions
    % This function is defined by N points on its half period (L/2)
    % This is equivalent to the Neumann problem on the full period (L)

    Lap = laplacian(1,N,'Neumann');
    D = derivative(1,N,'symmetric','Neumann');

    % kappa*(u_xx/u - u_x/u^2) + u - ln(u) + mu=0
    % dx = L/(2*N)

    bvp = @(u,mu,L) kappa*((2*N/L)^2*(Lap*u)./u - (2*N/L)*((D*u)./u).^2) ...
                        + u - log(u) + mu;
end