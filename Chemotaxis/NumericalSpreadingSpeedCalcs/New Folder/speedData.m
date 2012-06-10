chi = linspace(1.1,5);
for x=chi
    plot(x,speed(.5,x))
    drawnow
    hold on
end