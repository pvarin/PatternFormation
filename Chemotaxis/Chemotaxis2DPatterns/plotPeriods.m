close all

mu = -1.1;
kappa = 1;
b = 1;
for V0=linspace(.655,1.917)
    plot(V0,chemoPeriod(mu,kappa,b,V0),'o')
    hold on
    drawnow
end