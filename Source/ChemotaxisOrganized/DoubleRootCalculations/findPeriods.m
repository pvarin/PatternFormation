%% findPeriods.m
% Finds periods over a range of masses (massRange) for a given kappa using
% the double root computation for the Chemotaxis equations

function [masses L] = findPeriods(kappa, massRange)
    %% Define the dispersion relation
    % The dispersion relation is given by the linearized PDE:
    %       u_t = u_xx + s*u_x - mass*v_xx
    %       v_t = kappa*v_xx + s*v_x - v + u
    
    dr = @(nu,lambda,s,mass) ...
            (nu.^2+s.*nu-lambda).*...
            (kappa.*nu.^2+s.*nu-1-lambda) + ...
            mass.*nu.^2;
        
    dNu_dr = @(nu,lambda,s,mass) ...
            (2*nu+s).* (kappa.*nu.^2+s.*nu-1-lambda) + ...
            (nu.^2+s.*nu-lambda).*(2*kappa.*nu+s) + ...
            2*mass.*nu;
    
    % Define the problem that fsolve should solve
    problem = @(nu_r,nu_i,lambda_r,lambda_i,s,mass) ...
            [real(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             imag(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             real(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             imag(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass))];
    
    %% Begin with analytic solutions under the condition s=0
    % We are considering only the initial mass
    % We are looking for solutions where nu = i*k and lambda > 0
    % These can be found analytically
    
    s = 0;
    mass = massRange(1);
    
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
%     sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s,mass),sol)
    
    %% Follow the solution in s until Re(lambda) is close to 0
    
    opts = optimset('Display','off');
    
    sStep = 0.1;
    while sol(3)>0
        s = s+sStep;
        
        sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s,mass),sol,opts);
        fprintf('The real part of lambda is currently: %f\n',sol(3))
    end
   
    %% Follow the solution in mass and save the period data
    
    %restructure the solution vector
    sol = [sol(1) sol(2) sol(4) s];
    
    % initialize the mass and corresponding period vector
    masses = linspace(massRange(1), massRange(2))';
    L = zeros(length(masses),1);
    
    for i = 1:length(masses)
        sol = fsolve(@(a) problem(a(1),a(2),0,a(3),a(4),masses(i)),sol,opts);
        L(i) = 2*pi*sol(4)/sol(3);
    end
end