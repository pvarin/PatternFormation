function [V masses periods] = equilibriumSolutions(kappa,massRange)
    %% Constants and Parameters
    % Fundamental Constants
    global e
    e = exp(1);
    
    % Model Parameters (includes solving the double root problem to find
    % the correct period)
    mass = massRange(1);
    dblRootSol = solveChemDoubleRoot(kappa,mass);
    LDes = 2*pi*dblRootSol(4)/dblRootSol(3);
    
    % Secondary Parameters
    maxMu = maximumMu(kappa,mass,LDes);
    mu = maxMu-.1;
    v0 = 1;
    dx = 0.01;
    
    %% Create an Initial Guess
    % Find a solution to the differential equation with the initial
    % condition v=v0, v_x=0
    [v u] = full_chemIVPShoot(v0,dx,kappa,mu);
    N = length(v);
    L = dx*N*2;
    
    %% Create the Boundary Value Problem Solver
    bvp = chemBvp(kappa,N);
    
    chemSolveOpts = optimset('TolX',1e-8,'Display','off');
    chemSolve = @(v_,mu,L,bvp) fsolve(@(v) bvp(v,mu,L),v_,chemSolveOpts);
    
    %% Follow the solution in L in increments of ~1
    Lstep = (LDes-L)/floor((LDes-L)/1);
    
    while abs(LDes-L) > 1e-3
        L = L+Lstep;
        v = chemSolve(v,mu,L,bvp);
        displayState()
    end
    
    %% Follow the equilibrium solution in mu through to approximately the
    %% correct mass
    % Note: increasing mu will decrease the mass
    m = calcMass();
%     if m > 
    
    fprintf('Desired Mass: %f \n',mass)
    muStep = .1*sign(calcMass()-mass);%FIXME: the solutions seem to be fairly delicate in mu this step may need to be adapted
    
    while sign(calcMass()-mass)*muStep > 0
        displayState()
        if (mu + muStep > maxMu)
            warning('Mu is greater than the max!!!')
        end
        mu = mu + muStep;
        v = chemSolve(v,mu,L,bvp);
        plot(v);drawnow

    end
    fprintf('Final Mass: %f \n',calcMass())
    plot(v)
    
    %% Change the BVP to include the mass equation
    bvp = chemBvpMass(kappa,N);

    function [v mu] = chemSolveMass(v_guess,mu_guess,L,mass,bvp)
        opts = optimset('TolX',1e-8,'Display','off');
        guess = [v_guess; mu_guess];
        sol = fsolve(@(sol_) bvp(sol_,mass,L),guess,opts);
        v = sol(1:end-1);
        mu = sol(end);
    end

    %%  Follow the function in mass and log the progress
    fprintf('Solving something useful!! Do Work!!!\n')
    
    % Initialize the solution variables
    masses = linspace(massRange(1),massRange(2));
    periods = zeros(length(masses),1);
    V = zeros(N,length(masses));
    
    % Create the double root problem (to follow the appropriate period in
    % the mass)
    dbrProblem = doubleRootProblem(kappa);
    opts = optimset('Display','off');
    solveDbrProblem = @(guess,mass) fsolve(@(a) dbrProblem(a(1),a(2),0,a(3),a(4),mass),guess,opts);
    
    for i=1:length(masses)
        dblRootSol = solveDbrProblem(dblRootSol,masses(i));
        L = 2*pi*dblRootSol(4)/dblRootSol(3);
        [v mu] = chemSolveMass(v,mu,L,masses(i),bvp);
        
        periods(i) = L;
        V(:,i) = v;
        plot(linspace(L/N,L/2-L/N,N),v); drawnow
        displayState()
    end

    %% Nested Helper Functions
    function m=calcMass()
        m = mean(v);
    end

    function displayState()
        fprintf('\tPeriod: %6f Mass:%6f\n',L,calcMass())
    end
end