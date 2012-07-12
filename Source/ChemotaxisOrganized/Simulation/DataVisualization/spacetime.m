function spacetime(V,x,t)
    figure
    pcolor(x,t,V')
    shading flat
    colormap gray
    axis([0 x(end) 0 t(end) 0 .1 0 max(V(:,end/2))])
    drawnow
end