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
    
    s = 0;
    
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
%     sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s),sol)
    
    %% Follow the solution in s until Re(lambda) = 0
    
    opts = optimset('Display','off');
    
    % arrive close to the speed where the real part of lambda = 0
    sStep = 0.1;
    while sol(3)>0
        s = s+sStep;
        lastVal = sol(3);
        
        sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s),sol,opts);
        fprintf('The real part of lambda is currently: %f\n',sol(3))
    end
    
    % restructure the solution vector to give s as a parameter and solve
    % for real(lambda) = 0
    sol = [sol(1),sol(2),sol(4),s];
	sol = fsolve(@(a) problem(a(1),a(2),0,a(3),a(4)),sol,opts);
    
    %% Calculate the period
    % The period (L) is given by
    %   L = 2*pi*s/w
    % where w is the imaginary part of lambda
    
    L = 2*pi*sol(4)/sol(3);
    
    %% Notes
    %{
    Eventually we will follow this solution in mass, in order to relate the
    period and the mass when finding the periodic equilibria to the PDE
    
    ^This is implemented in findPeriods
    %}
end