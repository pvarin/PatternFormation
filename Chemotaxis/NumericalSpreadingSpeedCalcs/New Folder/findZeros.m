function y=findZeros(func)
    %fminsearch options
    options = optimset('fminsearch');
    options.TolX = 1e-6;
    
    %locate a good starting position
    x = fminsearch(@(a) -func(a,0),.01,options);
    y = fminsearch(@(b) func(x,b),.01,options);
    
    %iterative minimum approximation
    while (func(x,y)>0.1)
        x = fminsearch(@(a) func(a,y),x,options);
        y = fminsearch(@(b) func(x,b),y,options);
    end
    
    %linear interpolation iterative minimum approximation
    line1 = zeros(2); line2 = zeros(2);
    while (func(x,y)>1e-3)
        for i=1:2
            x = fminsearch(@(a) func(a,y),x,options);
            line1(i,:) = [x y];
            y = fminsearch(@(b) func(x,b),y,options);
            line2(i,:) = [x y];
        end
        coeff1 = polyfit(line1(end-1:end,1),line1(end-1:end,2),1);
        coeff2 = polyfit(line2(end-1:end,1),line2(end-1:end,2),1);
        x = (coeff1(2)-coeff2(2))/(coeff2(1)-coeff1(1));
        y = coeff1(1)*x+coeff1(2);
    end
end