load('Instability2Periods.mat')

plot(all_eVals)

unstableU = -all_eVects(1:end/2,2,100);
unstableV = -all_eVects(end/2+1:end,2,100);

stableU = all_eVects(1:end/2,3,100);
stableV = all_eVects(end/2+1:end,3,100);

subplot(2,1,1)
plot([unstableU unstableV])
subplot(2,1,2)
plot([stableU stableV])

figure
L = periods(100);
x = linspace(0,L,length(stableU))';
gamma = 2*pi/L;
rotateStableU = stableU.*exp(1i*gamma*x);
res = real(rotateStableU) - unstableU;

plot(res)
res'*res