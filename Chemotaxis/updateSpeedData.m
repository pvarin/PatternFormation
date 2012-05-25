function updateSpeedData(U,x,t,kappa,b)
    load('speedData.mat');
    speedData = [speedData; [spreadingSpeed(U,x,t) kappa b]];
    save('speedData.mat','speedData')
end