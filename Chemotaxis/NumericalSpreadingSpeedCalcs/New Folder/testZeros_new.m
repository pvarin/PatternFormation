clear all; close all;
kappa = 1;
chi = 1.01;

n = 40; m = 40;
real = linspace(0,.0005,n);
imag = linspace(0,.00005,m);
[s,L] = meshgrid(real,imag*1i);
wedgeProd = zeros(m,n);
for i = 1:n
    for j = 1:m
     wedgeProd(j,i) = eigenvectors(L(j,i),s(j,i),kappa,chi);
    end
end


close all
surf(real,imag,wedgeProd);
xlabel('s (real)');
ylabel('L (imaginary)');
axis square
drawnow
shg
hold on
func = @(a,b) eigenvectors(a*1i,b,kappa,chi);
findZeros_new_new(func)
