function sol = solveChemDoubleRoot(kappa,mass)
    
    problem = doubleRootProblem(kappa);

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
    
    %% Follow the solution in s until Re(lambda) ~ 0
    
    opts = optimset('Display','off');
    
    % arrive close to the speed where the real part of lambda = 0
    sStep = 0.1;
    sol = [0 k lambda 0];
    while sol(3)>0
        s = s+sStep;
        sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s,mass),sol,opts);
    end
    
    % restructure the solution vector to give s as a parameter and solve
    % for real(lambda) = 0
    sol = [sol(1),sol(2),sol(4),s];
	sol = fsolve(@(a) problem(a(1),a(2),0,a(3),a(4),mass),sol,opts);
end