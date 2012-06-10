clear all; close all;
n = 100; m = 100;
real = linspace(-5,5,n);
imag = linspace(-5,5,m);
[s,si] = meshgrid(real,imag*1i);
L = 3*1i;
s = s+si;
wedgeProd = zeros(m,n);
for i = 1:n
    for j = 1:m
     wedgeProd(j,i) = eigenvectors(L,s(j,i),.5,1.1);
    end
end
surf(real,imag,wedgeProd);
xlabel('real');
ylabel('imaginary');