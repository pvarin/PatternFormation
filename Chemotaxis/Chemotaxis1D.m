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

function varargout = Chemotaxis1D(varargin)
    %% Defaults
    %default model parameters
    b = 2;
    kappa = 1;

    %default simulation parameters
    dt = .01; numsteps = 10000;
    dx = .1; N = 1000;
    
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
            case 'SAVEWAVELENGTH'
                saveWavelength = varargin{2*i};
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
    
    nSaveStep = floor(numsteps/100);
    numsteps = nSaveStep*100;
    
    %% Prepare the Matrices
    % initialize the solution matrices
    U = zeros(N,1); V = zeros(N,1);
    U_ = zeros(N,floor(numsteps/nSaveStep)); V_ = zeros(N,floor(numsteps/nSaveStep));
    U = U0; V = V0;
    
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
        nstep = nstep+1;
        U = U_U\(L_U\(U_op*U-dt*b*D*(D*V.*U)));
        V = U_V\(L_V\(V_op*V+dt*U));
        %split the data into multiple parts
        if mod(nstep,nSaveStep)==0
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
        x = dx:dx:dx*N;
        t = nSaveStep*dt:nSaveStep*dt:dt*numsteps;
    end
        
    % Display Results
    if graph
        spacetime(U,x,t)   
        animate(U,V,dx)
    end
    
    % Save Data
    if saveData
        save('Data/CrankNicolsonChemotaxisData_new.mat','U','V','x','t','dx','dt','b','kappa')
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
            ['speed:      ', num2str(spreadingSpeed(U(5:end,:),x(5:end),t))]
            ['wavelength: ', num2str(wavelength(U(:,end),x))]
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