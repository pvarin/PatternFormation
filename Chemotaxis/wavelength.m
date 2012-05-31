function wavelength = wavelength(U, x)
    xi=[];
    dx = x(2)-x(1);
    mU = max(U);
    for i=2:length(U)-1
        if U(i)>max(U)/3 && U(i)>U(i+1) && U(i)>U(i-1)
            xi = [xi; i];
        end
    end
    wavelength = mean(diff(xi(2:end)*dx));
    
%     ['Found ' num2str(length(xi)) ' peaks']
%     xi*dx
end