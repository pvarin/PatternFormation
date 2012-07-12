function animate(C, E, dx)
    % unpack structure size
    N = size(C);
    numsteps = N(2);
    N = N(1);

    figure
    hold all
    for i=1:floor(numsteps/100):numsteps
        clf
        shg
        plot(linspace(0,dx*N,N),C(:,i))
        hold all
        plot(linspace(0,dx*N,N),E(:,i))
        axis([0 dx*N -2 2])
        drawnow
    end
end