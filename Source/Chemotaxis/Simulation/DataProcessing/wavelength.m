%% wavelength.m
%
% Given a vectror for U and a uniform mesh for x, this function calculates
% the wavelength of the pattern in U

function wavelength = wavelength(U, x)
    % initialize xi
    xi=[];
    
    % calculate the grid size and the maximum value of U
    dx = x(2)-x(1);
    mU = max(U);
    
    % search for local maxima that are at least 1/3 of the maximum value of
    % U
    for i=2:length(U)-1
        if U(i)>max(U)/3 && U(i)>U(i+1) && U(i)>U(i-1)
            xi = [xi; i];
        end
    end
    
    % drop the first peak (this one is likely adversely affected by
    % boundary conditions)
    wavelength = mean(diff(xi(2:end)*dx));
end