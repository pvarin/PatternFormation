load('../Data/equilibriumSolutions_kappa_10_mass_11_20.mat')

%two periods
evals = zeros(length(U(1,:)),4);

for i = 1:length(U(1,:));
    u2 = [U(:,i); U(end:-1:1,i); U(:,i); U(end:-1:1,i)];
    v2 = [V(:,i); V(end:-1:1,i); V(:,i); V(end:-1:1,i)];
    Lin = chemLinear(u2,v2,periods(i)*4,kappa);
    [~,eVals] = eigs(Lin,20,5);
    eVals = sort(diag(eVals),1,'descend');
    evals(i,:) = eVals(1:4);
    i
end

masses_ = masses;
load('../Data/equilibriumSolutions_kappa_10_mass_20_50.mat')
evals = [evals;zeros(length(U(1,:)),4)];

for i = 1:length(U(1,:));
    u2 = [U(:,i); U(end:-1:1,i); U(:,i); U(end:-1:1,i)];
    v2 = [V(:,i); V(end:-1:1,i); V(:,i); V(end:-1:1,i)];
    Lin = chemLinear(u2,v2,periods(i)*4,kappa);
    [~,eVals] = eigs(Lin,20,5);
    eVals = sort(diag(eVals),1,'descend');
    evals(i+length(masses_),:) = eVals(1:4);
    i+length(masses_)
end

masses = [masses_;masses];
plot(masses,evals','.')
xlabel('mass')
ylabel('eigenvalues')

save('Coarsening4Periods.mat','masses','evals')