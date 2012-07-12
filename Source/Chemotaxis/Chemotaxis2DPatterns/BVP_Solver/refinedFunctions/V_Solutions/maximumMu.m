%% maximumMu.m
% Only certain values of mu admit periodic solutions to the equilibrium
% chemotaxis equations. This function determines the maximum value of mu
% that admits periodic solutions

function mu = maximumMu(kappa,mass,T)
    temp = kappa.*(2*pi./T).^2+1;
    mu = min(log(mass)-mass,log(temp)-temp);
end