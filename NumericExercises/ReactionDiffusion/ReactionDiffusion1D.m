%% ReactionDiffusion1D.m
% Simulates reaction diffusion equations with specified initial conditions
% Neumann boundary conditions
% 
% ReactionDiffusion(IC, X)
%   IC: An NX2 matrix containing the initial conditions for the two
%       chemical species
%       e.g. [linspace(0,1)',linspace(2,4)']
% 
%   X:  A 2-vector containing the spatial boundaries points
%       e.g. [x_0, x_L]

function ReactionDiffusion1D(IC,X,a)
    
    %Extract initial conditions
    C0 = IC(:,1);
    E0 = IC(:,2);

    %Set simulation variables
    N = length(C0);
    nplotstep = 1;
    dt = .00001;
    dx = (X(2)-X(1))/N;
    kappa = .0001;
    gamma = 1;
    
    %Create the implicit Euler step matrix for Neumann boundary conditions
    uDiag = [0; ones(N-1,1)];
    cDiag = [-1;-2*ones(N-2,1);-1];
    lDiag = [ones(N-1,1); 0];
    B = [uDiag, cDiag, lDiag];
    d = [1; 0; -1];
    L = spdiags(B, d, N, N);
    
    S_C = speye(N) - dt/dx^2*L;
    S_C = S_C\speye(N);
    S_E = speye(N) - kappa*dt/dx^2*L;
    S_E = S_E\speye(N);
    [Top,Restart,Pause,Quit,CPlot,EPlot] = initgraphics(N,X);

    %run until the user exits
    while get(Quit,'value') == 0
        
        %initial conditions
        nstep = 0;
        C = C0;
        E = E0;
        
        set(Restart,'value',0)
        
        %loop until user interrupts
        while get(Restart,'value')==0 && get(Quit,'value')==0
            while get(Pause,'value')==1 && get(Quit,'value')==0
                pause(.1)
            end
            nstep = nstep+1;
            C = S_C*(C-F(C,E,gamma,a));
            E = S_E*(E+F(C,E,gamma,a));
            
            %plot the function at appropriate intervals
            if mod(nstep,nplotstep) == 0
                t = nstep*dt;
                set(CPlot,'ydata',C);
                set(EPlot,'ydata',E);
                set(Top,'string',sprintf('t = %6.2f',t))
                drawnow
            end
        end
    end
    close gcf
end

function [Top,Restart,Pause,Quit,CPlot,EPlot] = initgraphics(N,X)
   clf
   shg
   set(gcf,'numbertitle','off','name','ReactionDiffusion1D')
   subplot(2,1,1)
   CPlot = plot(linspace(X(1),X(2),N),zeros(1,N));
   axis([X(1) X(2) -.5 .5])
   Top = title('1-D Reaction Diffusion Simulation');
   subplot(2,1,2)
   EPlot = plot(linspace(X(1),X(2),N),zeros(1,N));
   axis([X(1) X(2) -2 2])
   Restart = uicontrol('position',[20 20 80 20],'style','toggle','string','restart');
   Pause = uicontrol('position',[120 20 80 20],'style','toggle','string','pause');
   Quit = uicontrol('position',[220 20 80 20],'style','toggle','string','stop');
end

function transfer = F(C,E,gamma,a)
    transfer = E.*(1-E).*(E-a) + gamma*C;
end