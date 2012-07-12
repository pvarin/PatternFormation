function [x_new F_new] = refineMeshArcLength(x,xi,xf,F,ds)
    N = length(F);

    % Assumes Homogeneous Neumann Boundary Conditions
    Lap = nonUniformLaplacian(x,xi,xf,'');
    D = nonUniformDerivative(x,xi,xf,'');
    
    F_x = D*F;
    F_xx = Lap*F;
    
    f = @(a) interp1(x,F,a);
    f_x = @(a) interp1(x,F_x,a);
    f_xx = @(a) interp1(x,F_xx,a);
    
    subplot(3,1,1);
    plot(x,f(x)); drawnow
    subplot(3,1,2);
    plot(x,f_x(x)); drawnow
    subplot(3,1,3);
    plot(x,f_xx(x)); drawnow
    
    x_new(1) = x(1);
    i=1;
    while x_new(end) < x(end)
%         dx = min(ds/sqrt(1+f_x(x_new(i))^2),ds_x/sqrt(1+f_xx(x_new(i))^2));
        dx = ds/sqrt(1+f_x(x_new(i))^2+f_xx(x_new(i))^2);
        x_new(i+1)=x_new(i)+dx;
        i=i+1;
    end
    
    F_new = f(x_new);
    
    figure% visualize the newly resampled function
    plot(x_new,F_new,'o'); drawnow
    figure% plot the mesh density
    plot(1./diff(x_new))
end