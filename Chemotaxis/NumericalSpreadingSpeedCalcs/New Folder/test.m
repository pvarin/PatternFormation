clear all; close all;
n = 50; m = 50;
real = linspace(-.1,.1,n);
imag = linspace(-.2,.2,m);
[L,Li] = meshgrid(real,imag*1i);
L = L+Li;
s = 0.2;
wedgeProd = zeros(m,n);
for i = 1:n
    for j = 1:m
     wedgeProd(j,i) = eigenvectors(L(j,i),s,.5,1.1);
    end
end
surf(real,imag,wedgeProd);
xlabel('real');
ylabel('imaginary');
% shading interp
% K = 1; X = 1;
% g= 4*(2*K^2*L.*(1+L)+K*L.*sqrt(K^3*L.*(1+L))+sqrt(K^3*L.*(1+L)).*(1+L-X))/K;
% plot(L,g,'-g')