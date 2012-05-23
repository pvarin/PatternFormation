x = (linspace(0,1).*sin(linspace(0,4*pi)))';
% HeatEqn1D(x,[0,10])
% ModifiedHeatEqn1D(x,[0,10],@F);

a=.5;
c = zeros(101,1);
e = ones(101,1)*a;
e(floor(end/2-3):floor(end/2+3)) = e(floor(end/2-3):floor(end/2+3))+.1;
ReactionDiffusion1D([c,e],[0,10],a)