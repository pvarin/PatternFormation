%% Chemotaxis1D.m
% Simulates and displays the voltage field on a 2 dimensional LC grid
% 
% TwoDimensionalGrid()
%   Runs the simulation with Neumann boundary conditions and arbitrary
%   predeteremined initial conditions on an arbitrary predeterimined
%   domain
%
% TwoDimensionalGrid(IC, X)
%   Runs the simulation with Neumann boundary conditions and an
%   initial condition specified by IC on the domain specified by X
%
%   IC: a 2XN matrix that specifies the initial conditions of the bacteria
%       concentration profile (first row) and the chemical concentration
%       profile (second row) of the system
%       e.g. [V0; U0]
%
%   X:  a 2 vector containing the lower and upper bounds of the spatial
%       domain
%       e.g. [0; 1]
%
% TwoDimensionalGrid(IC, X, identifier1, value1, identifier2, value2)
%   Runs the simulation with an initial condition specified by V0 and
%   options specified by identifier, value pairs
%
%   IC:         a matrix specifying the initial concentration profiles at
%               each point on the line
%   X:          a vector containing the bounds of the spatial domain
%   identifier: a string that specifies the option that will be set with
%               the proceeding value parameter
%   value:      the value of the option to be set
%   **a table of (identifier/value) pairs is provided below
%   
%   identifier          value
%
%   'Boundaries'        A string specifying the behavior of the boundary
%                       conditions
%                       'Neumann':  The simulation will enforce homogeneous
%                           Neumann boundary conditions Ux=Vx=0
%                       'Dirichlet':The simulation will enforce homogeneous
%                           Neumann boundary conditions U=V=0
%                       Default: 'Neumann'

function Chemotaxis1D(varargin)
    
    %set initial condition
    if isempty(varargin)
        %sample condition (for no inputs)
        X = linspace(0,1);
        V0 = sin(X).*X;
        U0 = cos(X).*X;
        X = [X(1),X(2)];
    else
        %extract initial condition
        IC = varargin{1};
        V0 = IC(1,:);
        U0 = IC(2,:);
        X = varargin{2};
    end
    
    %determine the size of the grid
    N=length(V0);%number of rows (column size)
    
    %set default values
    NEUMANN = 0;
    DIRICHLET = 1;
    boundaries = NEUMANN;
    
    %read variable inputs
    for i=1:floor((length(varargin)-1)/2)
        if strcmp(varargin{i*2},'Boundaries')
            if strcmp(varargin{i*2+1},'Neumann')
                boundaries = NEUMANN;
            elseif strcmp(varargin{i*2+1},'Dirichlet')
                boundaries = DIRICHLET;
            end
        end
    end
    
    %set simulation variables
    nplotstep = 10;
    dt = .1;
    [Top,Restart,Pause,Quit,VPlot,UPlot] = initgraphics(N,X);

    %run until the user exits
    while get(Quit,'value') == 0
        
        %initial conditions
        nstep = 0;
        V = [0 V0 0];
        U = [0 U0 0];
        dVdt = zeros(1,N+2);
        dUdt = zeros(1,N+2);
        set(Restart,'value',0)
        
        %loop until user interrupts
        while get(Restart,'value')==0 && get(Quit,'value')==0
            while get(Pause,'value')==1 && get(Quit,'value')==0
                pause(.1)
            end
            nstep = nstep+1;
            V = V + dVdt*dt;
            U = U + dUdt*dt;
            
            if boundaries == NEUMANN
                V(1) = V(2); U(1) = U(2);
                V(N) = V(N-1); U(N) = U(N-1);
            elseif boundaries == DIRICHLET
                V(1) = 0; U(1) = 0;
                V(N) = 0; U(N) = 0;
            end
            
            n=2:(N-1);
            dVdt(n) = (V(n-1)+V(n+1)-2*V(n));%FIXME: implement this
            dUdt(n) = (U(n-1)+U(n+1)-2*U(n));%FIXME: implement this
            
            %plot the function at appropriate intervals
            if mod(nstep,nplotstep) == 0
                t = nstep*dt;
                set(VPlot,'ydata',V(2:end-1));
                set(UPlot,'ydata',U(2:end-1));
                set(Top,'string',sprintf('t = %6.2f',t))
                drawnow
            end
        end
    end
    close gcf
end

function [Top,Restart,Pause,Quit,VPlot,UPlot] = initgraphics(N,X)
   clf
   shg
   set(gcf,'numbertitle','off','name','1-D Chemotaxis')
   subplot(2,1,1)
   VPlot = plot(linspace(X(1),X(2),N),zeros(1,N));
   axis([X(1) X(2) -2 2])
   Top = title('1-D Chemotaxis Simulation');
   subplot(2,1,2)
   UPlot = plot(linspace(X(1),X(2),N),zeros(1,N));
   axis([X(1) X(2) -2 2])
   Restart = uicontrol('position',[20 20 80 20],'style','toggle','string','restart');
   Pause = uicontrol('position',[120 20 80 20],'style','toggle','string','pause');
   Quit = uicontrol('position',[220 20 80 20],'style','toggle','string','stop');
end