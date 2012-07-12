load('../Data/equilibriumSolutions_kappa_10_mass_20_50.mat')

% Mass:  2
% Kappa: 1

u = [U(:,1); U(end:-1:1,1)];
v = [V(:,1); V(end:-1:1,1)];

Lin = chemLinear(u,v,periods(1)*2,kappa);

[eVect,eVal] = eigs(Lin,10,3);
eVal = diag(eVal)
[a, maxIndex] = max(eVal);
a

eigU = eVect(1:end/2,maxIndex);
eigV = eVect(end/2+1:end,maxIndex);
plot([eigU,eigV])
drawnow

% Mass:  5
% Kappa: 1

u = [U(:,100); U(end:-1:1,100)];
v = [V(:,100); V(end:-1:1,100)];

Lin = chemLinear(u,v,periods(100)*2,kappa);

[eVect,eVal] = eigs(Lin,10,3);
eVal = diag(eVal)
[a, maxIndex] = max(eVal);
a

eigU = eVect(1:end/2,maxIndex);
eigV = eVect(end/2+1:end,maxIndex);
figure
plot([eigU,eigV])
drawnow