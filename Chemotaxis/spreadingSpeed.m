function speed = spreadingSpeed(U, x, t)
    [~, Y]  = gradient(U');
        l = max(Y(:))/2;
    s = size(Y);
    xi = []; yi=[];
    for j=1:s(2)
        [m i] = max(Y(:,j));
        if m>l
            xi = [xi; i];
            yi = [yi; j];
        end
    end
    p = polyfit(yi,xi,1);
    drawnow
    speed = p(1)*(x(2)-x(1))/(t(2)-t(1));

    pcolor(Y); shading interp
    hold on
    plot(yi,xi,'*g')
    x_ = linspace(0,s(1));
    plot(x_,p(2)+p(1)*x_,'r')
end