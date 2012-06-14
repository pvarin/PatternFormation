kappa = 4;
chi = 2;

a1 = @(k) ((kappa+1)*k.^2+1)/2;
a2 = @(k) ((kappa-1)*k.^2+1)/2;
L = @(k) -a1(k)+sqrt(a2(k).^2+chi*k.^2);
dLdk = @(k) (-kappa-1)*k+(a2(k)*(kappa-1).*k+chi*k)./sqrt(a2(k).^2+chi*k.^2);

k = linspace(-1,1);
plot(k,L(k))
hold all
% plot(k(1:end-1),diff(L(k))./diff(k))
plot(k,zeros(size(k)),'k')
plot(k,dLdk(k))

% plot(k,((kappa+1)^2*(a2(k).^2+chi*k.^2)-a2(k).^2*(kappa-1)^2-2*a2(k)*chi*(kappa-1)-chi^2)/100)
% plot(k,(a2(k).^2*4*kappa-2*a2(k)*chi*(kappa-1)+(kappa+1)^2*chi*k.^2-chi^2)/100)

a = kappa*(kappa-1)^2;
b = 4*chi*kappa+2*kappa*(kappa-1);
c = kappa-chi*(kappa-1)-chi^2;

plot(k,(a*k.^4+b*k.^2+c)/100)

k_max1 = sqrt((-b+sqrt(b^2-4*a*c))/2);
k_max2 = (-b-sqrt(b^2-4*a*c))/2;

% plot(k_max1,L(k_max1),'o')
% plot(k_max2,L(k_max2),'o')