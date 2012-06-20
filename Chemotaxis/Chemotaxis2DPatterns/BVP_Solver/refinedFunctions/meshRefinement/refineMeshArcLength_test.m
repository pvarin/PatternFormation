function refineMeshArcLength_test(kappa,V,L,ds,ds_x)
    N = length(V);
    dx = L/(2*N);
    x = 0:dx:(N-1)*dx;
    x = x + dx/2;
    
    E = ones(N,1);
    D = spdiags([-E E],[-1 1],N,N);
    D(1,1) = -1; D(end,end) = 1;
    D = D/(2*dx);
    V_x = D*V;

    Lap = spdiags([E -2*E E],-1:1,N,N);
    Lap(1,1) = -1; Lap(end,end) = -1;
    Lap = Lap/dx^2;
    V_xx = Lap*V;
    
    
    f = @(a) interp1(x,V,a);
    f_x = @(a) interp1(x,V_x,a);
    f_xx = @(a) interp1(x,V_xx,a);
    
    plot(x,f(x)); drawnow; hold all
    plot(x,f_x(x)); drawnow
    plot(x,f_xx(x)); drawnow
    
    x_new(1) = dx/2;
    i=1;
    while x_new(end) < x(end)
        dx = min(ds/sqrt(1+f_x(x_new(i))^2),ds_x/sqrt(1+f_xx(x_new(i))^2));
        x_new(i+1)=x_new(i)+dx;
        i=i+1;
    end
    figure
    plot(x_new,f(x_new),'o'); drawnow
    figure
    plot(1./diff(x_new))
    length(x_new)
end