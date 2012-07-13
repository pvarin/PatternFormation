%% animate2D.m
%
% Takes a file path to a file containing variables 'U', 'x', 'y', and 't', 
% where U contains spacetime data and x, y, and t specify the mesh for this
% spacetime data

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