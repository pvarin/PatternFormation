load('../Data/equilibriumSolutions_kappa_10_mass_20_50.mat')

% Mass:  2
% Kappa: 1

k = linspace(0,5,100);
% k = [0 1]';
% k = 1;
% k = 0;

% one period
numVals1 = 4;
evals1 = zeros(length(k),numVals1);

u = U(:,1);
v = V(:,1);

Lin1 = chemLinear(u,v,periods(1),kappa);
T1 = transverseTerms(u,v);

for i=1:length(k)
    Lin_transverse = Lin1 + k(i).^2.*T1;
    [~,eVals] = eigs(Lin_transverse,10,6);
    eVals = sort(diag(eVals),1,'descend');
    evals1(i,:) = eVals(1:numVals1);
    i
end

figure(1)
plot(k,evals1,'o')

% two periods
numVals2 = 6;
evals2 = zeros(length(k),numVals2);

u = [u; u(end:-1:1)];
v = [v; v(end:-1:1)];

Lin2 = chemLinear(u,v,periods(1)*2,kappa);
T2 = transverseTerms(u,v);

for i=1:length(k)
    Lin_transverse = Lin2 + k(i)^2*T2;
    [~,eVals] = eigs(Lin_transverse,10,6);
    eVals = sort(diag(eVals),1,'descend');
    evals2(i,:) = eVals(1:numVals2);
    i
end

figure(2)
plot(k,evals2,'o')