% x = (logspace(0,1,500)-1)/9*pi;
% Li = x(1);
% Lf = x(end);
% x = x(2:end-1)';

x = .1:.2:99.9;
x = x';
Li = 0;
Lf = 100;

Deriv = nonUniformDerivative(x,Li,Lf,'Dirichlet');

[V D] = eig(full(Deriv));
[V D] = sortem(V,D);
plot(x,V(:,end-4:end))
D = diag(D);
D(end-4:end)