%% spreadingSpeed.m
%
% Given a spacetime diagram for U as well as a uniform mesh for x and t,
% this function calculates the spreading speed of patterns in U

function speed = spreadingSpeed(U, x, t)
    % calculate the time derivative of U
    [~, Y]  = gradient(U');
    
    % find the maximum time derivative on the spacetime diagram
    l = max(Y(:))/3;
    
    s = size(Y);
    ti = []; xi=[];
    dx = x(2)-x(1);
    dt = t(2)-t(1);
    
    for j=1:s(2)
        % find the maximum time derivative for each point in time if this
        % derivative is greater then 1/3 of the absolute maximum, then
        % store the time and space indices
        [m i] = max(Y(:,j));
        if m>l
            ti = [ti; i];
            xi = [xi; j];
        end
    end
    
    % fit a line to the peaks
    p = polyfit(ti*dt,xi*dx,1);
    % the speed is the slope of this line
    speed = p(1);
    
%     figure
%     pcolor(x,t,Y); shading interp
%     hold on
%     plot(xi*dx,ti*dt,'*g')
%     x_ = linspace(0,s(1)*dx);
%     plot(x_,(x_-p(2))/speed,'r')
end