load('../../Data/equilibriumSolutions_kappa_10_mass_11_50.mat')

% Mass:  5
% Kappa: 1

% one period
u = U(:,200);
v = V(:,200);
L = periods(200);

k = 2.5;

Lin = chemLinear(u,v,L,kappa);
T = transverseTerms(u,v);

Lin_transverse = Lin + k^2*T;
[eVect, eVals] = eigs(Lin_transverse,10,6);
[eVect, eVals] = sortEigen(eVect,eVals);

eVals(1)
plot([eVect(1:end/2,1) eVect(end/2+1:end,1)]);

% two periods
u = [u; u(end:-1:1)];
v = [v; v(end:-1:1)];
L = L*2;

Lin = chemLinear(u,v,L,kappa);
T = transverseTerms(u,v);

Lin_transverse = Lin + k^2*T;
[eVect, eVals] = eigs(Lin_transverse,10,6);
[eVect, eVals] = sortEigen(eVect,eVals);

eVals(1)
eVals(2)
figure
plot([eVect(1:end/2,1) eVect(end/2+1:end,1)]);
figure
plot([eVect(1:end/2,2) eVect(end/2+1:end,2)]);
