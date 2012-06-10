clear all; close all;
kappa = .5;
chi = 2;

n = 100; m = 100;
real = linspace(0,5,n);
imag = linspace(0,5,m);
[s,L] = meshgrid(real,imag*1i);
wedgeProd = zeros(m,n);
for i = 1:n
    for j = 1:m
     wedgeProd(j,i) = eigenvectors(L(j,i),s(j,i),kappa,chi);
    end
end
%%
surf(real,imag,wedgeProd);
xlabel('s (real)');
ylabel('L (imaginary)');
axis equal
drawnow

hold on

options = optimset('fminsearch');
options.TolX = 1e-6;

func = @(a,b) eigenvectors(a*1i,b,kappa,chi);
x = fminsearch(@(a) -func(a,0),.01,options);
y = fminsearch(@(b) func(x,b),.01,options);
while (func(x,y)>0.1)
    x = fminsearch(@(a) func(a,y),x,options);
    plot3(y,x,func(x,y),'r*'); drawnow
    y = fminsearch(@(b) func(x,b),y,options);
    plot3(y,x,func(x,y),'r*'); drawnow
end
line1 = zeros(2); line2 = zeros(2);
while (func(x,y)>1e-3)
    for i=1:2
        x = fminsearch(@(a) func(a,y),x,options);
        plot3(y,x,func(x,y),'g*'); drawnow
        line1(i,:) = [x y];
        y = fminsearch(@(b) func(x,b),y,options);
        line2(i,:) = [x y];
        plot3(y,x,func(x,y),'g*'); drawnow
    end
    coeff1 = polyfit(line1(:,1),line1(:,2),1);
    coeff2 = polyfit(line2(:,1),line2(:,2),1);
    
%     x = (coeff1(2)-coeff2(2))/(coeff2(1)-coeff1(1));
%     y = coeff1(1)*x+coeff1(2);

    x1 = fminsearch(@(a) func(a,coeff1(1)*a+coeff1(2)), x,options);
    x2 = fminsearch(@(a) func(a,coeff1(1)*a+coeff1(2)), x,options);
    x = mean([x1,x2]);
    y1 = coeff1(1)*x1+coeff1(2);
    y2 = coeff2(1)*x2+coeff2(2);
    y = mean([y1,y2]);
    plot3(y,x,0,'b*'); drawnow
end
func(x,y)