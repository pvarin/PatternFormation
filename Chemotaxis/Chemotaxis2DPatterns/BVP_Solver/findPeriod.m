function L = findPeriod(kappa, mass)
    %% Define the dispersion relation
    % The dispersion relation is given by the linearized PDE:
    %       u_t = u_xx + s*u_x - mass*v_xx
    %       v_t = kappa*v_xx + s*v_x - v + u
    
    dr = @(nu,lambda,s) ...
            (nu.^2+s.*nu-lambda).*...
            (kappa.*nu.^2+s.*nu-1-lambda) + ...
            mass.*nu.^2;
        
    dNu_dr = @(nu,lambda,s) ...
            (2*nu+s).* (kappa.*nu.^2+s.*nu-1-lambda) + ...
            (nu.^2+s.*nu-lambda).*(2*kappa.*nu+s) + ...
            2*mass.*nu;
    
    % Define the problem that fsolve should solve
    problem = @(nu_r,nu_i,lambda_r,lambda_i,s) ...
            [real(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s));
             imag(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s));
             real(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s));
             imag(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s))];
    
    %% Begin with analytic solutions under the condition s=0
    % We are looking for solutions where nu = i*k and lambda > 0
    % These can be found analytically
    
    % Solve for nu
    a = kappa.*(kappa-1).^2;
    b = (2*kappa.*(kappa-1)+4*mass.*kappa);
    c = kappa - mass.*(kappa-1) - mass.^2;
    k = sqrt((-b+sqrt(b.^2-4*a.*c))./(2*a));
    
    % Solve for lambda
    a1 = ((kappa+1)*k.^2+1)/2;
    a2 = ((kappa-1)*k.^2+1)/2;
    lambda = -a1+sqrt(a2.^2+mass*k.^2);
    
    %Confirm that this is a solution
    sol = [0 k lambda 0];
    sol = fsolve(@(nu_r,nu_i,l_r,l_i) problem(nu_r,nu_i,l_r,l_i,0),sol)
    
    %% Follow the solution in s until Re(lambda) = 0
    
    %% Notes
    %{
    Eventually we will follow this solution in mass, in order to relate the
    period and the mass when finding the periodic equilibria to the PDE
    %}
end