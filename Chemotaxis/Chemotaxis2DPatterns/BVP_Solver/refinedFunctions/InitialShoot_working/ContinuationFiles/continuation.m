function continuation(kappa, massRange)
    %% Setup
    masses = linspace(massRange(1), massRange(2));

    if ~exist('initialSolution.mat','file')
        warning('initialSolution.mat does not exist you should make one here')
        return
    end
    % initial solution should have the variables u, v, L, kappa, mu
    load('initialSolution.mat');
    
    %extract data
    N = length(u);
    dx = L/N;
    mass = mean(u);
    
    %determine the direction to take
    [~ massIndex] = abs(mass-masses);
    startingMassIndex = massIndex;
    
    %determine the selected period
    dblRootSol = solveChemDoubleRoot(kappa,masses(massIndex));
    LDes = 2*pi*dblRootSol(4)/dblRootSol(3);
    %% Follow the Solution in the Period
    
end