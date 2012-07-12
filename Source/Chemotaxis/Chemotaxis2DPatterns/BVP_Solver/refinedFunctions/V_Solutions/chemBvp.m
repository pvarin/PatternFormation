function bvp=chemBvp(kappa,N)
    % Defines boundary value problem as function in N points using finite differences
    % kappa v_xx+ exp(mu+v)-v=0

    % We are solving for even periodic solutions
    % This function is defined by N points on its half period (L/2)
    % This is equivalent to the Neumann problem on the full period (L)

    E=ones(N,1);
    Lap=spdiags([E -2*E E],-1:1,N,N);
    Lap(1,1)=-1; Lap(N,N)=-1;

    % k*V_xx = V - exp(u+b*V)
    % dx = L/(2*N)

    bvp = @(v,mu,L) kappa/(L/(2*N))^2*Lap*v + exp(mu+v) - v;
end