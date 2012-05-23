x = (linspace(0,1).*sin(linspace(0,4*pi)))';
% HeatEqn1D(x,[0,10])
% ModifiedHeatEqn1D(x,[0,10],@F);

a=.5;
c = zeros(100,1);
e = ones(100,1)*a;
e(floor(end/2)) = e(floor(end/2))+.01;
ReactionDiffusion1D([c,e],[0,1],a)