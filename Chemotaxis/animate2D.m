function animate2D(path)
    load(path)
    figure
    f = pcolor(x,y,U(:,:,1));
    shading interp
    colormap gray
    axis([0 x(end) 0 y(end) 0 3 0 2])
    for i=1:length(t)
        set(f,'cdata',U(:,:,i))
        drawnow
    end
end