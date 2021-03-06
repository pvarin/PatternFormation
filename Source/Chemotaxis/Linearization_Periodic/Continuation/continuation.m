function continuation(kappa, massRange)
    %% Setup
    masses = linspace(massRange(1), massRange(2))';
    initialDataFile = sprintf('../../../../Data/PeriodicInitialGuess/initialSolution_kappa_%d.mat',floor(10*kappa));
    if ~exist(initialDataFile,'file')
        %FIXME: implement this
        warning('initialSolution.mat does not exist you should make one here')
        return
    end
    
    % initial solution should have the variables u, v, L, (kappa? this
    % might conflict with the specified kappa), (mu?)
    load(initialDataFile);
    
    %extract data
    N = length(u);
    mass = mean(u);
    
    %create the boundary value problem
    bvp = chemBVP(N);
    
    %verify that this is a solution to the BVP
    sol = [u;v;0;0];
    
    plot(sol)
    figure
    plot(bvp(sol,L,kappa,mass))

    %determine the closest mass
    [~, massIndex] = min(abs(mass-masses));
    
    %determine the selected period
    dblRootSol = solveChemDoubleRoot(kappa,masses(massIndex));
    LDes = 2*pi*dblRootSol(4)/dblRootSol(3);
    

    
    %% Obtain the correct mass
    sol = [u;v;0;0];
    figure(1)
    plot([sol(1:N),sol(N+1:2*N)])
    drawnow
    shg
    fprintf('Obtaining the correct mass...\n')
    fprintf('\tDesired Mass: %f\n',masses(massIndex))
    fprintf('\tCurrentMass: %f\n',mass)
    mass = masses(massIndex);
    sol = chemSolve(sol,L,kappa,mass,bvp);
    figure(1)
    plot([sol(1:N),sol(N+1:2*N)])
    drawnow
    shg
    
    %% Follow the Solution in the Period
    fprintf('Following the solution in the period...\n')
    fprintf('\tDesired Period: %f\n',LDes)
    Lstep = (LDes-L)/floor((LDes-L)/.1);
    
    while abs(LDes-L) > 1e-3
        fprintf('\tCurrent Mass: %f\n\tCurrent Period: %f\n',mean(sol(1:N)),L);
        L = L+Lstep;
        fprintf('\tSolving for a period of %f...\n',L)
        sol = chemSolve(sol,L,kappa,mass,bvp);
        figure(2)
        plot([sol(1:N),sol(N+1:2*N)])
        drawnow
    end
    fprintf('Final Mass: %f\nFinal Period: %f\n',mean(sol(1:N)),L);
    
    save(initialDataFile,'kappa','u','v','L')
    
    %% Follow the Solution in the Mass
    fprintf('\nFollowing the solution in the mass...\n')
    U = zeros(N,length(masses));
    V = zeros(N,length(masses));
    periods = zeros(length(masses),1);
    
    U(:,massIndex) = sol(1:N);
    V(:,massIndex) = sol(N+1:2*N);
    periods(massIndex) = L;
    
    for i=massIndex-1:-1:1
        fprintf('\tComputing double root...\n')
        dblRootSol = solveChemDoubleRoot(kappa,masses(i));
        L = 2*pi*dblRootSol(4)/dblRootSol(3);
        
        %step in the period
%         fprintf('Solving for mass: %f period: %f\n',masses(i+1),L)
%         sol = chemSolve(sol,L,kappa,masses(i+1),bvp);
        %step in the mass
        fprintf('Solving for mass: %f period: %f\n',masses(i),L)
        sol = chemSolve(sol,L,kappa,masses(i),bvp);
        %store the data
        U(:,i) = sol(1:N);
        V(:,i) = sol(N+1:2*N);
        periods(i) = L;
        
        %plot the data to show progress
        figure(2)
        plot([sol(1:N),sol(N+1:2*N)])
        drawnow
        fprintf('\tMass Dummy: %f\n\tShift Dummy: %f\n',sol(end-1),sol(end))
    end
    
    sol(1:2*N) = [U(:,massIndex); V(:,massIndex)];
    
    for i=massIndex+1:length(masses)
        fprintf('\tComputing double root...\n')
        dblRootSol = solveChemDoubleRoot(kappa,masses(i));
        L = 2*pi*dblRootSol(4)/dblRootSol(3);
        
        %step in the mass and the period simultaneously
        fprintf('Solving for mass: %f period: %f\n',masses(i),L)
        sol = chemSolve(sol,L,kappa,masses(i),bvp);
        %store the data
        U(:,i) = sol(1:N);
        V(:,i) = sol(N+1:2*N);
        periods(i) = L;
        
        %plot the data to show progress
        figure(2)
        plot([sol(1:N),sol(N+1:2*N)])
        drawnow
        fprintf('\tMass Dummy: %f\n\tShift Dummy: %f\n',sol(end-1),sol(end))
    end
    
    filename = sprintf('../../../../Data/PeriodicEquilibria/equilibriumSolutions_kappa_%d_mass_%d_%d.mat',floor(10*kappa),floor(10*masses(1)),floor(10*masses(end)));
    save(filename,'U','V','kappa','masses','periods')
end