%% chemLinear.m
%
% Returns the matrix of the linearization of the chemotaxis equations, Lin,
% about the equilibrium u0, v0
%
% The linearized chemotaxis equations are then given by:
%   U_t = Lin*U
%
% Requires the equilibrium solutions u0 and v0 as column vectors as well as
% the size of the domain, L, and the parameter kappa

function Lin = chemLinear(u0,v0,L,kappa)
    %% Construct the Derivative and Laplacian Matrices
    N = length(u0);
    dx = L/N;
    D = derivative(dx,N,'symmetric','periodic');
    Lap = laplacian(dx,N,'periodic');
    
    %% Compute Derivatives
    v0_x = D*v0;
    v0_xx = D*v0_x;
    
    u0_x = D*u0;
    
    %% Construct the Linearization in Parts
    u = spdiags(u0,0,N,N);
    u_x = spdiags(u0_x,0,N,N);
    
    v_x = spdiags(v0_x,0,N,N);
    v_xx = spdiags(v0_xx,0,N,N);
    
    I = speye(N);
    
    Lin = [Lap-v_xx-v_x*D,   -u_x*D-u*D*D;
           I,                kappa*Lap-I];
end