function V = equilibriumSolutions(kappa,massRange)
    %% Constants and Parameters
    % Fundamental Constants
    global e
    e = exp(1);
    
    % Model Parameters
    mass = massRange(1);
    dblRootSol = initDoubleRootSolution(kappa,mass);
    LDes = 2*pi*dblRootSol(4)/dblRootSol(3)
    
    % Secondary Parameters
    maxMu = maximumMu(kappa,mass,LDes);
    mu = maxMu-.1;
    v0 = 1;
    dx = 0.01;
    
    %% Create an Initial Guess
    % Find a solution to the differential equation
    
    %% Nested Helper Functions
    function [v N L]=initialShoot(v0,dx)
        % Solves the IVP 
        % Used to generate an inital guess for V
        v(1)=v0;
        v(2)=v(1);
        j=2;
        while (v(j) >=v(j-1))
            v(j+1)=2*v(j)-v(j-1)+dx^2/kappa*(-exp(mu+v(j))+v(j));
            j=j+1;
        end
        N=j-1;
        L=N*dx*2;
        v=v(2:length(v))';
    end

    function vsol=chemSolve(initial,bvp)
        % Finds the solutions to the boundary value problem defined in bvp
        % using fsolve

        options = optimset('TolX',1e-8,'Display','off');
        vsol=fsolve(@(v) bvp(v,mu,L),initial,options);
    end

    function m=calcMass()
        m = mean(v);
    end    
end

function mu = maximumMu(kappa,mass,T)
    temp = kappa.*(2*pi./T).^2+1;
    mu = min(log(mass)-mass,log(temp)-temp);
end

function sol = initDoubleRootSolution(kappa,mass)
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
    
    %Confirm that this is a solution
    sol = [0 k lambda 0];
%     sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s),sol)
    
    %% Follow the solution in s until Re(lambda) = 0
    
    opts = optimset('Display','off');
    
    % arrive close to the speed where the real part of lambda = 0
    sStep = 0.1;
    
    while sol(3)>0
        s = s+sStep;
        sol = fsolve(@(a) problem(a(1),a(2),a(3),a(4),s,mass),sol,opts);
    end
    
    % restructure the solution vector to give s as a parameter and solve
    % for real(lambda) = 0
    sol = [sol(1),sol(2),sol(4),s];
	sol = fsolve(@(a) problem(a(1),a(2),0,a(3),a(4),mass),sol,opts);
end

function drp = doubleRootProblem(kappa)
    dr = @(nu,lambda,s,mass) ...
            (nu.^2+s.*nu-lambda).*...
            (kappa.*nu.^2+s.*nu-1-lambda) + ...
            mass.*nu.^2;
        
    dNu_dr = @(nu,lambda,s,mass) ...
            (2*nu+s).* (kappa.*nu.^2+s.*nu-1-lambda) + ...
            (nu.^2+s.*nu-lambda).*(2*kappa.*nu+s) + ...
            2*mass.*nu;
    
    % Define the problem that fsolve should solve
    drp = @(nu_r,nu_i,lambda_r,lambda_i,s,mass) ...
            [real(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             imag(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             real(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             imag(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass))];

end