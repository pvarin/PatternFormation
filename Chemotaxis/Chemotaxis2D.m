%% Chemotaxis2D.m
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
%   'U0'        -The initial condition for U, this must be an MxN matrix
%               that corresponds to the initial concentration profile of
%               the bacteria
%               default: [1.1; ones(N-1,1)]
%   'V0'        -The initial condition for V, this must be an MxN matrix
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
%   'dx'        -The grid cell size in the x direction, i.e. the distance 
%               between consecutive horizontal grid points, this must be a 
%               positive scalar
%               default: .1
%   'dy'        -The grid cell size in the y direction, i.e. the distance 
%               between consecutive vertical grid points, this must be a
%               positive scalar
%               default: .1
%   'numsteps'  -The number of timesteps to take, this must be a positive
%               integer
%               default: 10000
%   'N'         -The grid size in the x direction, i.e. the number of
%               horizontal grid cells, this must be a positive integer
%               default: 1000
%   'M'         -The grid size in the y direction, i.e. the number of
%               vertical grid cells, this must be a positive integer
%               default: 1000
%   'graph'     -Option to animate the simulation results
%               default: true
%   'saveData'  -Option to save the results of the simulation in a file
%               named 'CrankNicolsonChemotaxisData_new.mat', this must be a
%               Boolean value
%               default: false

function varargout = Chemotaxis2D(varargin)
    %% Defaults
    %default model parameters
    b = 2;
    kappa = 1;

    %default simulation parameters
    dt = .01; numsteps = 10000;
    dx = .1; N = 1000;
    dy = .1; M = 1000;
    
    U0init = false; V0init = false;
    Ninit = false; Minit = false;
    saveData = false;
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
                U0s = size(U0);
                if Ninit && N ~= U0s(2)
                    err = MException('Input:Overdefined', ...
                    'The value of N that you specified is inconsistent with the dimensions of U0');
                    throw(err)
                end
                if Minit && M ~= U0s(1)
                    err = MException('Input:Overdefined', ...
                    'The value of M that you specified is inconsistent with the dimensions of U0');
                    throw(err)
                end
                N = U0s(2);
                M = V0s(1);
                clear V0s U0s
            case 'V0'
                V0 = varargin{2*i};
                V0init = true;
                if U0init && min(size(U0) ~= size(V0))
                    err = MException('Input:Overdefined', ...
                    'The dimensions of U0 and V0 must be consistent');
                    throw(err)
                end
                V0s = size(U0);
                if Ninit && N ~= V0s(2)
                    err = MException('Input:Overdefined', ...
                    'The value of N that you specified is inconsistent with the dimensions of U0');
                    throw(err)
                end
                if Minit && M ~= V0s(1)
                    err = MException('Input:Overdefined', ...
                    'The value of M that you specified is inconsistent with the dimensions of U0');
                    throw(err)
                end
                N = U0s(2);
                M = V0s(1);
                clear V0s U0s
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
            case 'DY'
                dy = varargin{2*i};
            case 'NUMSTEPS'
                numsteps = varargin{2*i};
            case 'N'
                N = varargin{2*i};
                Ninit = true;
                if U0init
                U0s = size(U0);
                end
                if V0init
                    V0s = size(V0);
                end
                if (U0init && U0s(2) ~= N) || (V0init && V0s(2) ~= N)
                    err = MException('Input:Overdefined', ...
                    'The value of N that you specified is inconsistent with the dimensions of the initial conditions');
                    throw(err)
                end
                clear U0s V0s
            case 'M'
                M = varargin{2*i};
                Minit = true;
                if U0init
                    U0s = size(U0);
                end
                if V0init
                    V0s = size(V0);
                end
                if (U0init && U0s(1) ~= M) || (V0init && V0s(1) ~= M)
                    err = MException('Input:Overdefined', ...
                    'The value of M that you specified is inconsistent with the dimensions of the initial conditions');
                    throw(err)
                end
                clear U0s V0s
            case 'SAVEDATA'
                saveData = varargin{2*i};
            case 'GRAPH'
                graph = varargin{2*i};
            otherwise
                warning(['Argument ' num2str(2*i-1) ' is not a valid parameter name'])
        end
    end
    
    %default initial conditions
    if ~U0init
        U0 = ones(M,N);
        U0(floor(end/2)-2:floor(end/2)+2,floor(end/2)-2:floor(end/2)+2) = U0(floor(end/2)-2:floor(end/2)+2,floor(end/2)-2:floor(end/2)+2) + .2;
    end
    if ~V0init
        V0 = ones(M,N);
    end
    
    nSaveStep = floor(numsteps/100);
    numsteps = nSaveStep*100;
    
    if graph
        [Top, UPlot] = initgraphics(N,M,dx,dy);
    end
    
    %% Prepare the Matrices
    % initialize the solution matrices
    U = U0; V = V0;
    
    if saveData
        U_ = zeros(M,N,floor(numsteps/nSaveStep)); V_ = zeros(M,N,floor(numsteps/nSaveStep));
    end
    
    
    %prepare the operator matrices
    Dx = derivativeMatrix(N,dx); Dy = derivativeMatrix(M,dy);
    Lx = laplacianMatrix(N,dx); Ly = laplacianMatrix(M,dy);

    [L_Diff_Ux, U_Diff_Ux, ~] = lu(speye(N)-dt/2*Lx);
    [L_Diff_Uy, U_Diff_Uy, ~] = lu(speye(M)-dt/2*Ly);
    [L_Diff_Vx, U_Diff_Vx, ~] = lu(speye(N)-dt/2*kappa*Lx);
    [L_Diff_Vy, U_Diff_Vy, ~] = lu(speye(M)-dt/2*kappa*Ly);

    Diffuse_Ux = speye(N)+dt/2*Lx;
    Diffuse_Uy = speye(M)+dt/2*Ly;
    Diffuse_Vx = speye(N)+dt/2*kappa*Lx;
    Diffuse_Vy = speye(M)+dt/2*kappa*Ly;


    %% Simulate
    %Crank-Nicolson Method
    nstep = 1;
    while nstep<=numsteps
        if graph
            set(UPlot,'zdata',U)
            set(Top,'string',sprintf('t = %6.2f',dt*nstep))
            drawnow
        end
        nstep = nstep+1;
        %Crank-Nicolson Diffusion
        U = ((U*Diffuse_Ux)/U_Diff_Ux)/L_Diff_Ux;
        V = ((V*Diffuse_Vx)/U_Diff_Vx)/L_Diff_Vx;
          
        U = U_Diff_Uy\(L_Diff_Uy\(Diffuse_Uy*U));
        V = U_Diff_Vy\(L_Diff_Vy\(Diffuse_Vy*V));
        
        %Crank-Nicolon Other Linear
        V = V*(1+dt/2)/(1-dt/2);
        
        %Explicit Euler Non Linear and Cross Terms
        U = U + dt*b*(V*Dx'.*U)*Dx';
        U = U + dt*b*Dy*(Dy*V.*U);
        
        V = V + U*dt;
%         U = ((U*Ux_op-dt*b*(V*Dx.*U)*Dx)/U_Ux)/L_Ux;
%         V = ((V*Vx_op+dt*U)/U_Vx)/L_Vx;
        
%         U = U_Uy\(L_Uy\(Uy_op*U-dt*b*Dy*(Dy*V.*U)));
%         V = U_Vy\(L_Vy\(Vy_op*V+dt*U));
        %split the data into multiple parts
        if saveData && mod(nstep,nSaveStep)==0
            i=nstep/nSaveStep;
            U_(:,:,i)=U;
            V_(:,:,i)=V;
        enddrawnow
        end
    end
    
    %% Save Data
    if saveData
        U = U_;
        V = V_;
        clear U_ V_
        
        x = dx:dx:dx*N;
        y = dy:dy:dy*M;
        t = nSaveStep*dt:nSaveStep*dt:dt*numsteps;
        
        save('Data/2DChemotaxisData_new.mat','U','V','x','y','t','dx','dy','dt','b','kappa')
    end
    
    %% Construct output
    switch nargout
        case 0
            %do nothing
        case 2
            varargout = {U, V};
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

function [Top,UPlot] = initgraphics(N,M,dx,dy)
   clf
   shg
   set(gcf,'numbertitle','off','name','2-D Chemotaxis Simulation')
   UPlot = surf(dx:dx:N*dx,dy:dy:M*dy,zeros(M,N));
   shading interp
   axis([0 dy*M 0 dx*N 0 5])
   xlabel('x')
   ylabel('y')
   Top = title('2-D Chemotaxis Simulation');
end