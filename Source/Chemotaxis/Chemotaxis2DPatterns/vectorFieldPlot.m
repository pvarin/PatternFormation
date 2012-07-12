mu = -2;

[V Vx] = meshgrid(linspace(0,5,20), linspace(-5,5,21));
[U Ux] = meshgrid(linspace(.6,5,20), linspace(-5,5,21));

Vxx = V - exp(mu+V);
U = exp(V+mu);
Ux = Vx.*U;
Uxx = (Vxx+Vx.^2).*U;

quiver(V,Vx,Vx,Vxx)
figure
quiver(U,Ux,Ux,Uxx)