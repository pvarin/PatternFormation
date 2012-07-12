%% HeatEqn1D.m
% Simulates the heat equation with a specified initial condition and
% Neumann boundary conditions using an implicit Euler's method
% 
% HeatEqn1D(U0, X)
%   U0: A column vector of N elements containing the initial condition for
%       the heat equation
%       e.g. linspace(0,1)'
%   X:  A 2-vector containing the spatial boundaries points
%       e.g. [x_0, x_L]

function HeatEqn1D(U0,X)
    
    %set simulation variables
    N = length(U0);
    nplotstep = 10;
    dt = .001;
    dx = (X(1)-X(2))/N;
    L = circshift(eye(N),[0,1])+ circshift(eye(N),[0,-1]) + -2*eye(N);
    %Neumann Boundary Conditions
    L(1,1) = -1; L(end,end) = -1;
    L(1,end) = 0; L(end,1) = 0;
    S = eye(N) - dt/dx^2*L;
    S = S^-1;
    
    [Top,Restart,Pause,Quit,UPlot] = initgraphics(N,X);

    %run until the user exits
    while get(Quit,'value') == 0
        
        %initial conditions
        nstep = 0;
        U = U0;
        set(Restart,'value',0)
        
        %loop until user interrupts
        while get(Restart,'value')==0 && get(Quit,'value')==0
            while get(Pause,'value')==1 && get(Quit,'value')==0
                pause(.1)
            end
            nstep = nstep+1;
            U = S*U;
            
            %plot the function at appropriate intervals
            if mod(nstep,nplotstep) == 0
                t = nstep*dt;
                set(UPlot,'ydata',U);
                set(Top,'string',sprintf('t = %6.2f',t))
                drawnow
            end
        end
    end
    close gcf
end

function [Top,Restart,Pause,Quit,UPlot] = initgraphics(N,X)
   clf
   shg
   set(gcf,'numbertitle','off','name','1-D Heat Equation')
   UPlot = plot(linspace(X(1),X(2),N),zeros(1,N));
   axis([X(1) X(2) -2 2])
   Top = title('1-D Heat Equation Simulation');
   Restart = uicontrol('position',[20 20 80 20],'style','toggle','string','restart');
   Pause = uicontrol('position',[120 20 80 20],'style','toggle','string','pause');
   Quit = uicontrol('position',[220 20 80 20],'style','toggle','string','stop');
end