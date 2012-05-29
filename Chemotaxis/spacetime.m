function spacetime(U,x,t)
    figure
    pcolor(x,t,U')
    shading flat
    colormap gray
    axis([0 x(end) 0 t(end) 0 .1 0 2])
    drawnow
end