x = (logspace(0,1,500)-1)/9*pi;
Li = x(1);
Lf = x(end);
x = x(2:end-1)';

Lap = nonUniformLaplacian(x,Li,Lf,'');

[V D] = eigs(Lap,4,'sm');
plot(x,V)
diag(D)