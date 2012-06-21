%% full_ChemBvp.m
%
% Defines the boundary value problem that should be solved later with
% fsolve.
%
% Specifically, this function returns a function handle that fsolve can
% use. The function that is returns takes in three parameters and returns
% the error:
%   sol: the solution, a column vector of the form [U; V; sigma]
%        at an exact solution this function will return the zero vector
%
%   L:   the length of the domain, this is necessary to recalculate dx each
%        time the bvp is called

function bvp = full_ChemBvp(kappa,N)
    
    Lap = laplacian(1,N,'Neumann');
    % The derivative left derivative with no boundary conditions applied to
    % the right derivative with Neumann boundary conditions is equivalent
    % to the Laplacian with Neumann boundary conditions
    %   To test try: 
    %     full(derivative(1,5,'left','')*derivative(1,5,'right','Neumann')) 
    D_l = derivative(1,N,'left','');
    D_r = derivative(1,N,'right','Neumann');
    
    % Extract U, V and calculate dx
    U = @(sol) sol(1:(end-1)/2);
    V = @(sol) sol((end+1)/2:end-1);
    dx = @(L,N) L/(2*N);
    
    % Create the boundary value problem
    bvp = @(sol,mass,L) ...
        [Lap/dx(L,N)^2*U(sol) - ...
            D_l/dx(L,N)*(U(sol).*(D_r/dx(L,N)*V(sol))) + ...
            sol(end);
            
         kappa*Lap/dx(L,N)^2*V(sol) + ...
            U(sol) - V(sol);
            
         mean(U(sol)) - mass];
end