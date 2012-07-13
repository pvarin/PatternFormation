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
eVal = diag(eVal);
[a, maxIndex] = max(eVal);

eigU = eVect(1:end/2,maxIndex);
eigV = eVect(end/2+1:end,maxIndex);
plot([eigU,eigV])