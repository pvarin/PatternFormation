function speed = spreadingSpeed(U, x, t)
    [~, Y]  = gradient(U');
    l = max(Y(:))/3;
    s = size(Y);
    ti = []; xi=[];
    dx = x(2)-x(1);
    dt = t(2)-t(1);
    for j=1:s(2)
        [m i] = max(Y(:,j));
        if m>l
            ti = [ti; i];
            xi = [xi; j];
        end
    end
    p = polyfit(ti*dt,xi*dx,1);
    drawnow
    speed = p(1);
    
%     figure
%     pcolor(x,t,Y); shading interp
%     hold on
%     plot(xi*dx,ti*dt,'*g')
%     x_ = linspace(0,s(1)*dx);
%     plot(x_,(x_-p(2))/speed,'r')
end