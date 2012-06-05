clear all; close all;
n = 20; m = 20;
real = linspace(-20,20,n);
imag = linspace(-10,10,m);
[L,Li] = meshgrid(real,imag*1i);
L = L+Li;
s = 5;
wedgeProd = zeros(m,n);
for i = 1:n
    for j = 1:m
     wedgeProd(j,i) = eigenvectors(L(j,i),s,.5,1.1);
    end
end
surf(real,imag,wedgeProd);
xlabel('real');
ylabel('imaginary');
% K = 1; X = 1;
% g= 4*(2*K^2*L.*(1+L)+K*L.*sqrt(K^3*L.*(1+L))+sqrt(K^3*L.*(1+L)).*(1+L-X))/K;
% plot(L,g,'-g')