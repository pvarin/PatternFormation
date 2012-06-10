function x = newsolve(func,guess)

    f_tol = 0.0001;
    d_tol = [0.0001;0.0001];
    x = guess;
    f_x = func(x);
    gradient = NaN;
    while(f_x>f_tol)
        if ~isnan(gradient)
            d_tol = (-f_x.*gradient)/100;
        end
        %calculate the gradient
        f_x_up = func(x+d_tol);
        gradient = (f_x_up-f_x)./d_tol;
        %update guess and function values
        if f_x>func(-f_x./gradient/2+x)
        x = -f_x./gradient/2+x
        f_x = func(x);
    end
end