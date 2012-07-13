%% ChemotaxisAnalysisGraphs.m
%
% Used to explore the instability of the 1-D Chemotaxis equations.
% When this script is run it will plot the dispersion relation of the
% linearized system for a variety of b (i.e. chi or average mass) values
% from this graph we can see that the system is not unstable until b>1 at
% which point only a finite range of wavenumbers is unstable, suggesting
% that the full non-linear system will exhibit some sort of wavenumber
% selection.

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