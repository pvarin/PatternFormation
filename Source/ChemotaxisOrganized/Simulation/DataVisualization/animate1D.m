function animate1D(U, V, dx)
    % unpack structure size
    N = size(U);
    numsteps = N(2);
    N = N(1);

    figure
    hold all
    for i=1:floor(numsteps/100):numsteps
        clf
        shg
        plot(linspace(0,dx*N,N),U(:,i))
        hold all
        plot(linspace(0,dx*N,N),V(:,i))
        axis([0 dx*N 0 max(U(:,end))])
        drawnow
    end
end