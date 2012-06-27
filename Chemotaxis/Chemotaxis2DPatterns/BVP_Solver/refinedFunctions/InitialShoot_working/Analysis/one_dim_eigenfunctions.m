load('../Data/equilibriumSolutions_kappa_10_mass_20_50.mat')

% Mass:  2
% Kappa: 1

k = 2;

u = U(:,1);
v = V(:,1);

Lin1 = chemLinear(u,v,periods(1),kappa);
T1 = transverseTerms(u,v);

Lin_transverse = Lin1 + k^2*T1;

[eVect,eVal] = eigs(Lin_transverse,10,5);
diag(eVal)