%% updateWavelengthData(U,x,kappa,b)
%
% Updates the data stored in the file Data/wavelengthData.mat

function updateWavelengthData(U,x,kappa,b)
    load('Data/wavelengthData.mat');
    wavelengthData = [wavelengthData; [wavelength(U(:,end),x) kappa b]];
    save('Data/wavelengthData.mat','wavelengthData')
end