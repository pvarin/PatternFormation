%% simulateChemotaxis1D.m
% Simulates the one dimensional chemotaxis equations:
%   U_t = U_xx - b(UV_x)_x
%   V_t = kappa*V_xx + U - V
%
% Returns the spacetime data for the bacteria concentration (U) and the
% chemical concentration (V)
%
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

function varargout = simulateChemotaxis1D(varargin)
    %import dataVisualization directory
    addpath('./DataVisualization')
    
    %% Defaults
    %default model parameters
    
    N = 500;
    kappa = 1;
     
    %default initial conditions
    U0 = 5*ones(N,1);
    V0 = 5*ones(N,1);
    
    %default simulation parameters
    dt = .01; numsteps = 3200;
    N = length(U0);
    dx = .1;
    
    %perturb the equilibrium
    epsU = zeros(N,1);
    epsV = [.1*ones(1,1); zeros(N-1,1)];
    
    U0 = U0 + epsU;
    V0 = V0 + epsV;
    
    U0init = false; V0init = false; Ninit = false;
    saveData = false; saveSpeed = false; saveWavelength = false;
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
                dataFilename = varargin{2*i};
            case 'SAVESPEED'
                saveSpeed = varargin{2*i};
            case 'SAVEWAVELENGTH'
                saveWavelength = varargin{2*i};
            case 'GRAPH'
                graph = varargin{2*i};
            otherwise
                warning(['Argument ' num2str(2*i-1) ' is not a valid parameter name'])
        end
    end
    
    if numsteps>1000
        nSaveStep = floor(numsteps/1000);
        numsteps = nSaveStep*1000;
    else
        nSaveStep = 1;
    end
    
    %% Prepare the Matrices
    % initialize the solution matrices
    U = zeros(N,1); V = zeros(N,1);
    U_ = zeros(N,floor(numsteps/nSaveStep)); V_ = zeros(N,floor(numsteps/nSaveStep));
    U = U0; V = V0;
    
    %prepare the operator matrices
    D = derivative(dx,N,'symmetric','Neumann');
    Lap = laplacian(dx,N,'Neumann');

    [L_U, U_U, ~] = lu(speye(N)-dt*Lap);
    [L_V, U_V, ~] = lu(speye(N)-dt*(kappa*Lap));

%     U_op = speye(N)+dt/2*Lap;
%     V_op = speye(N)+dt/2*(kappa*Lap-speye(N));

    %% Simulate
    %Semi-Implicit Method
    nstep = 1;
    while nstep<=numsteps
        nstep = nstep+1;
        U = U_U\(L_U\(U - dt*(D*((D*V).*U))));
        V = U_V\(L_V\(V + dt*(U - V)));
        %split the data into multiple parts
        if mod(nstep,nSaveStep)==0
            figure(1)
            plot([U,V])
            nstep
            drawnow
            i=nstep/nSaveStep;
            U_(:,i)=U;
            V_(:,i)=V;
        end
    end
    
    U = U_;
    V = V_;
    clear U_ V_
    
    %% Display Results and Save Data
    if graph || saveData || saveSpeed || saveWavelength
        x = dx/2:dx:dx*N;
        i = floor(linspace(1,floor(length(x)/2),2000));
        i = floor(linspace(1,length(x),2000));
        x = x(i);
        U = U(i,:);
        V = V(i,:);
        t = nSaveStep*dt:nSaveStep*dt:dt*numsteps;
    end
        
    % Display Results
    if graph
        spacetime(V,x,t)   
        animate1D(U,V,dx)
    end
    
    % Save Data
    if saveData
        U = cast(U,'single');
        V = cast(V,'single');
        x = cast(x,'single');
        t = cast(t,'single');
        save(dataFilename,'U','V','x','t','dx','dt','kappa')
    end
    if saveSpeed
        updateSpeedData(U,x,t,kappa,b)
    end
    if saveWavelength
        updateWavelengthData(U,x,kappa,b)
    end
    
    %% Construct output
    switch nargout
        case 0
            % do nothing
        case 1
            s = spreadingSpeed(U(5:end,:),x(5:end),t);
            varargout = {s};
        case 2
            varargout = {U, V};
        case 3
            s = spreadingSpeed(U(5:end,:),x(5:end),t);
            varargout = {s, U, V};
        otherwise
            err = MException('Output:Unsupported', ...
                    ['This function is unsupported for ' ...
                    num2str(nargout) ' output arguments']);
            throw(err)
    end
end