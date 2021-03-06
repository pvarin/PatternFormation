function refineMesh
    x_ = linspace(0,2*pi);
    v = sin(x_);
    
    ds = 0.1;
    
    g = max(x_);
    l = min(x_);
    f = @(a) interp1(x_,v,min(g,max(a,l)));
    f_x = @(a) (f(a+.001)-f(a-.001))./(min(a+.001,g)-max(a-.001,l));
    
    plot(x_,f(x_)); drawnow; hold all
    plot(x_,f_x(x_)); drawnow
    
    x(1) = 0;
    i=1;
    while x(end) < x_(end)
        x(i+1)=x(i)+ds/sqrt(1+f_x(x(i))^2);
        i=i+1;
    end
    hold all
    plot(x,f(x),'o'); drawnow
    figure
    plot(x)
end