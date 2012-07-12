%% updateWavelengthData(U,x,kappa,b)
%
% Updates the data stored in the file Data/speedData.mat

function updateSpeedData(U,x,t,kappa,b)
    load('Data/speedData.mat');
    speedData = [speedData; [spreadingSpeed(U,x,t) kappa b]];
    save('Data/speedData.mat','speedData')
end