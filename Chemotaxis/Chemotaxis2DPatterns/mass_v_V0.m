kappa = 1;
b = 1;
mu = -2;

V0 = linspace(0.355,4.18,100);
mass = zeros(size(V0));
period = zeros(size(V0));
for i = 1:length(V0)
    [period(i) mass(i)] = IVP_onePeriod(kappa,b,mu,V0(i));
end

plot(V0,mass,'.')
figure
plot(V0,period,'.')
figure
plot(period,mass,'.')