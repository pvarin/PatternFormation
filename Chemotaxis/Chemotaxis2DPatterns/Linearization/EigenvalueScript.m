load('EquilibriumSolutions_kappa_0-5.mat')

v = V(:,10);
mass = M(10);
period = L(10);

[err A] = linearizedOperator(0.5,v,period,'Periodic',2);

% b = eig(full(A));
% b = sort(real(b),1,'descend');
% b(1:6)

[eigVect eigVal] = eigs(A,5,-.1);
eigVal
err
close all
j=2;
subplot(2,1,1)
plot(eigVect(1:end/2,j))
subplot(2,1,2)
plot(eigVect(end/2+1:end,j))