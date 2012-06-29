fileSuffix = {'_mass_5-0000',...
              '_mass_3-4848',...
              '_mass_2-8788',...
              '_mass_2-7273',...
              '_mass_2-0000'};
          
dataPrefix = 'CoarseningChemotaxisData';
dataExtension = 'mat';
figurePrefix = '../Figures/CoarseningFigure';
figureExtension = 'fig';
     
for i=1:length(fileSuffix)
    load([dataPrefix fileSuffix{i} '.' dataExtension])
    figure(i)
    pcolor(x,t,cast(V','double'))
    shading interp
    colormap gray
    axis([0 x(end) 0 t(end) 0 0.00001 0 max(V(:,1))])
    drawnow
%     saveas(gcf, [figurePrefix fileSuffix{i}], figureExtension )
end