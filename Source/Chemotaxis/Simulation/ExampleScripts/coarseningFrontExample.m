addpath('..')
%% Aggregate Coarsening
%%Set Parameters
% model parameters
mass = 20;
kappa = 10;
% simulation parameters
N = 2760*2;
dx = .05;
numsteps =10000;
dt = .01;

%%Set Initial Condition
% homogeneous equilibrium
U0 = mass*ones(N,1);
V0 = U0;
% perturb the equilibrium
V0(1:10,1) = V0(1:10,1) + .1;

%%Run the Simulation
simulateChemotaxis1D('kappa',kappa,'U0',U0,'V0',V0,...
                     'N',N,'dx',dx,'numsteps',numsteps,'dt',dt)
                 
%% Parasitic Coarsening
%%Set Parameters
% model parameters
mass = 2;
kappa = 1;
% simulation parameters
N = 2760;
dx = .1;
numsteps =300000;
dt = .01;

%%Set Initial Condition
% homogeneous equilibrium
U0 = mass*ones(N,1);
V0 = U0;
% perturb the equilibrium
V0(1:10,1) = V0(1:10,1) + 1;

%%Run the Simulation
simulateChemotaxis1D('kappa',kappa,'U0',U0,'V0',V0,...
                     'N',N,'dx',dx,'numsteps',numsteps,'dt',dt)