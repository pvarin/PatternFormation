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

function varargout = Chemotaxis1D(varargin)
    %% Defaults
    %default model parameters
    load('equilibriumSolutions_kappa_10_mass_11_50.mat',...
         'U','V','periods','kappa')
    
    %default initial conditions
    U0 = [U(:,150); U(end:-1:1,150)];
    V0 = [V(:,150); V(end:-1:1,150)];
    
    %default simulation parameters
    dt = .01; numsteps = 80000;
    N = length(U0); dx = periods(150)*2/N;
    
    %perturb the equilibrium
    dU0_dm = [U(:,102); U(end:-1:1,102)];
    dV0_dm = [V(:,102); V(end:-1:1,102)];
    dU0_dm = dU0_dm - U0; dU0_dm = dU0_dm/norm(dU0_dm);
    dV0_dm = dV0_dm - V0; dV0_dm = dV0_dm/norm(dV0_dm);
    
    dU0_dx = U0([2:N 1])-U0([N 1:N-1]);
    dV0_dx = V0([2:N 1])-V0([N 1:N-1]);
    dU0_dx = dU0_dx/norm(dU0_dx);
    dV0_dx = dV0_dx/norm(dV0_dx);
    
    epsU = (dU0_dm + dU0_dx).*cos(linspace(0,2*pi,N))';
    epsV = (dV0_dm + dV0_dx).*cos(linspace(0,2*pi,N))';
    U0 = U0 + epsU;
    V0 = V0 + epsV;
    
    %shift the IC so that the peaks will coalesce in the center
    U0 = U0([N/4+1:N,1:N/4]);
    V0 = V0([N/4+1:N,1:N/4]);
    
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
    D = derivative(dx,N,'symmetric','periodic');
    Lap = laplacian(dx,N,'periodic');

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
        t = nSaveStep*dt:nSaveStep*dt:dt*numsteps;
    end
        
    % Display Results
    if graph
        spacetime(U,x,t)   
        animate(U,V,dx)
    end
    
    % Save Data
    if saveData
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

%% Helper Functions
function animate(U, V, dx)
    % unpack structure size
    N = size(U);
    numsteps = N(2);
    N = N(1);

    figure
    hold all
    for i=1:floor(numsteps/100):numsteps
        clf
        shg
        plot(linspace(0,dx*N,N),U(:,i))
        hold all
        plot(linspace(0,dx*N,N),V(:,i))
        axis([0 dx*N 0 25])
        drawnow
    end
end

function spacetime(U,x,t)
    figure
    pcolor(x,t,U')
    shading flat
    colormap gray
    axis([0 x(end) 0 t(end) 0 .1 0 2])
    drawnow
end