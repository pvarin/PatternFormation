%% Chemotaxis1D.m
% Simulates the one dimensional chemotaxis equations:
%   U_t = U_xx - b(UV_x)_x
%   V_t = kappa*V_xx + U - V
%
% Returns the spacetime data for the bacteria concentration (U) and the
% chemical concentration (V)
% [U,V] = Chemotaxis1D('parameter1', value1, 'parameter2', value2, ...)
%   -This function takes in a variable number of parameters, if a parameter
%   is not specified, the simulation will run with a predetermined default
%   value
%
% Parameters:
%   'U0'        -The initial condition for U, this must be a column vector
%               that corresponds to the initial concentration profile of
%               the bacteria
%               default: [1.1; ones(N-1,1)]
%   'V0'        -The initial condition for V, this must be a column vector
%               that corresponds to the initial concentration profile of
%               the chemical
%               default: ones(N,1)
%   'b'         -The chemotactic sensitivity, this must be a positive
%               scalar
%               default: 1.01
%   'kappa'     -The relative diffusifity of the chemical in comparison
%               with to the bacteria, this must be a positive scalar
%               default: 1
%   'dt'        -The simulation timestep, i.e. time between consecutive
%               simulation steps, this must be a positive scalar
%               default: .01
%   'dx'        -The grid cell size, i.e. the distance between consecutive
%               grid cells, this must be a positive scalar
%               default: .1
%   'numsteps'  -The number of timesteps to take, this must be a positive
%               integer
%               default: 10000
%   'N'         -The grid size, i.e. the number of grid cells, this must be
%               a positive integer
%               default: 1000
%   'graph'     -Option to plot the spacetime data and animate the
%               simulation results
%               default: true
%   'saveData'  -Option to save the results of the simulation in a file
%               named 'CrankNicolsonChemotaxisData_new.mat', this must be a
%               Boolean value
%               default: false
%   'saveSpeed' -Option to add the speed data to the file 'speedData.mat',
%               this must be a Boolean value
%               default: false

function [U,V] = Chemotaxis1D(varargin)
    %% Defaults
    %default model parameters
    b = 1.01;
    kappa = 1;

    %default simulation parameters
    dt = .01; numsteps = 10000;
    dx = .1; N = 1000;
    nSaveStep = 1000;
    
    U0init = false; V0init = false; Ninit = false;
    saveData = false; saveSpeed = false;
    graph = true;

    %% Input and Error handling
    for i=1:floor(nargin/2)
        switch upper(varargin{2*i-1})
            % initial conditions
            case 'U0'
                U0 = varargin{2*i};
                U0init = true;
                if V0init && min(size(U0) ~= size(V0))
                    err = MException('Input:Overdefined', ...
                    'The dimensions of U0 and V0 must be consistent');
                    throw(err)
                end
                if Ninit && N ~= length(U0)
                    err = MException('Input:Overdefined', ...
                    'The value of N that you specified is inconsistent with the dimensions of U0');
                    throw(err)
                end
                N = length(V0);
            case 'V0'
                V0 = varargin{2*i};
                V0init = true;
                if U0init && min(size(U0) ~= size(V0))
                    err = MException('Input:Overdefined', ...
                    'The dimensions of U0 and V0 must be consistent');
                    throw(err)
                end
                if Ninit && N ~= length(V0)
                    err = MException('Input:Overdefined', ...
                    'The value of N that you specified is inconsistent with the dimensions of V0');
                    throw(err)
                end
                N = length(V0);
            % model parameters
            case 'B'
                b = varargin{2*i};
            case 'KAPPA'
                kappa = varargin{2*i};
            % simulation parameters
            case 'DT'
                dt = varargin{2*i};
            case 'DX'
                dx = varargin{2*i};
            case 'NUMSTEPS'
                numsteps = varargin{2*i};
            case 'N'
                N = varargin{2*i};
                Ninit = true;
                if (U0init && length(U0) ~= N) || (V0init && length(V0) ~= N)
                    err = MException('Input:Overdefined', ...
                    'The value of N that you specified is inconsistent with the dimensions of the initial conditions');
                    throw(err)
                end
            case 'SAVEDATA'
                saveData = varargin{2*i};
            case 'SAVESPEED'
                saveSpeed = varargin{2*i};
            case 'GRAPH'
                graph = varargin{2*i};
            otherwise
                warning(['Argument ' num2str(2*i-1) ' is not a valid parameter name'])
        end
    end
    
    %default initial conditions
    if ~U0init
        U0 = ones(N,1);
        U0(1:5) = U0(1:5) - .2;
    end
    if ~V0init
        V0 = ones(N,1);
    end
    
    %% Prepare the Matrices
    % initialize the solution matrices
    U = zeros(N,nSaveStep); V = zeros(N,nSaveStep);
    U_ = zeros(N,numsteps); V_ = zeros(N,numsteps);
    U(:,1) = U0; V(:,1) = V0;
    
    %prepare the operator matrices
    D = derivativeMatrix(N,dx);
    L = laplacianMatrix(N,dx);

    [L_U, U_U, ~] = lu(speye(N)-dt/2*L);
    [L_V, U_V, ~] = lu(speye(N)-dt/2*(kappa*L-speye(N)));

    U_op = speye(N)+dt/2*L;
    V_op = speye(N)+dt/2*(kappa*L-speye(N));

    %% Simulate
    %Crank-Nicolson Method
    nstep = 1;
    while nstep<=numsteps
        i = mod(nstep,nSaveStep)+1;
        nstep = nstep+1;
        if i>1
            U(:,i) = U_U\(L_U\(U_op*U(:,i-1)-dt*b*D*(D*V(:,i-1).*U(:,i-1))));
            V(:,i) = U_V\(L_V\(V_op*V(:,i-1)+dt*U(:,i-1)));
        elseif i==1
            U(:,i) = U_U\(L_U\(U_op*U(:,end)-dt*b*D*(D*V(:,end).*U(:,end))));
            V(:,i) = U_V\(L_V\(V_op*V(:,end)+dt*U(:,end)));
        end
        %split the data into multiple parts
        if i==nSaveStep
            U_(:,nSaveStep*(floor(nstep/i)-1)+1:nstep)=U;
            V_(:,nSaveStep*(floor(nstep/i)-1)+1:nstep)=V;
        end
    end
    
    U = U_;
    V = V_;
    clear U_ V_
    
    %% Lower Data Resolution
    x = dx:dx:dx*N;
    if (numsteps > 1000) && (graph || saveData || saveSpeed)
        Dt = floor(numsteps/1000);
        U = U(:,1:Dt:numsteps); V = V(:,1:Dt:numsteps);
        x = linspace(0,dx*N,N);
        t = Dt*dt:Dt*dt:numsteps*dt;
        spacetime(U,linspace(0,dx*N,N),Dt*dt:Dt*dt:numsteps*dt)
    else
        t = dt:dt:dt*numsteps;
    end
    
    %% Display Results
    if graph
        spacetime(U,x,t)   
        animate(U,V,dx)
    end

    ['speed: ', num2str(spreadingSpeed(U(5:end,:),x(5:end),t))]
    
    %% Save Data
    if saveData || saveSpeed
        if saveData
            save('CrankNicolsonChemotaxisData_new.mat','U','V','x','t','dx','dt','b','kappa')
        end
        if saveSpeed
            updateSpeedData(U,x,t,kappa,b)
        end
    end
end

%% Helper Functions
function D = derivativeMatrix(N, dx)
    B = [-ones(N,1) ones(N,1)];
    d = [-1; 1];
    D = spdiags(B, d, N, N);
    D(1,1) = -1;
    D(N,N) = 1;
    D = D/(2*dx);
end

function L = laplacianMatrix(N, dx)
    uDiag = [0; ones(N-1,1)];
    cDiag = [-1;-2*ones(N-2,1);-1];
    lDiag = [ones(N-1,1); 0];
    B = [uDiag, cDiag, lDiag];
    d = [1; 0; -1];
    L = spdiags(B, d, N, N)/dx^2;
end