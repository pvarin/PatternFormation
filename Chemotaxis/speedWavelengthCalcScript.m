%% Varying b
% b = 1.2:.1:10;
% numsteps = ceil(interp1([1.2,2.3,2.7,5,10],[100000, 9000, 6000, 3000, 1500],b));
% N = ceil(interp1([1.2,2.3,2.7,5,10],[1000, 802, 730, 500,400],b));
% 
% 
% for i=1:length(b)
%     b(i)
%     Chemotaxis1D('b',b(i),'kappa',.5,'numsteps',numsteps(i),'N',N(i),'saveSpeed',true,'graph',false);
% end
% load('Data/speedData.mat')
% close all
% plot(speedData(:,3),speedData(:,1),'.')

%% Varying kappa
% kappa = .15:.1:5;
% numsteps = ceil(interp1([.15,1,5],[8000,18000,60000],kappa));
% N = ceil(interp1([.15,1,5],[500,1000,2000],kappa));
% 
% 
% for i=1:length(b)
%     kappa(i)
%     Chemotaxis1D('b',2,'kappa',kappa(i),'numsteps',numsteps(i),'N',N(i),'saveSpeed',true,'graph',false);
% end
% load('Data/speedData.mat')XS
% close all
% plot(speedData(:,2),speedData(:,1),'.')





%% Varying b
% b = 1.2:.05:4.5;
% numsteps = ceil(interp1([1.2,1.3,2.3,2.7,5,10],[200000, 40000, 9000, 6000, 3000, 1500],b));
% N = ceil(interp1([1.2,1.3,2.3,2.7,5,10],[2000, 1000, 802, 730, 500,400],b));
% 
% 
% for i=1:length(b)
%     b(i)
%     Chemotaxis1D('b',b(i),'kappa',.5,'numsteps',numsteps(i),'N',N(i),'saveWavelength',true,'graph',false);
% end
% load('Data/wavelengthData.mat')
% close all
% plot(wavelengthData(:,3),wavelengthData(:,1),'.')

%% Varying kappa
kappa = .15:.1:3.5;
numsteps = ceil(interp1([.15,1,5],[8000,18000,60000],kappa));
N = ceil(interp1([.15,1,5],[500,1000,2000],kappa));


for i=1:length(kappa)
    kappa(i)
    Chemotaxis1D('b',2,'kappa',kappa(i),'numsteps',numsteps(i),'N',N(i),'saveWavelength',true,'graph',false);
end
load('Data/wavelengthData.mat')
close all
plot(wavelengthData(:,2),wavelengthData(:,1),'.')
