k=linspace(-4,4);
kappa = 2;
b = 2;
u0 = 1;

param_vec = linspace(0,10);
[k b] = meshgrid(k,param_vec);

lambda_1 = (-(1+k.^2+kappa.*k.^2)+sqrt((1+k.^2+kappa.*k.^2).^2 - 4*(kappa.*k.^4+k.^2-b.*u0.*k.^2)))./2;
lambda_2 = (-(1+k.^2+kappa.*k.^2)-sqrt((1+k.^2+kappa.*k.^2).^2 - 4*(kappa.*k.^4+k.^2-b.*u0.*k.^2)))./2;

figure(1)
surf(k,param_vec, lambda_1);
xlabel('k')
ylabel('b')
% figure(2)
% surf(k,param_vec, lambda_2);

% figure(3)
% plot(k,[lambda_1; lambda_2])