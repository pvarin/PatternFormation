%% full_chemIVPShoot_symmetric.m
%
% Solves the initial value problem for the equilibrium states of the
% chemotaxis equations:
%
% 0 = U_xx - (U*V_x)_x
% 0 = V_xx - V + U
%
% There are strict condition on the parameters in order to obtain periodic
% orbits:
%   kappa > 0
%   mu < -1
%   v0 must also be chosen wisely to obtain a periodic orbit
%       *note that v0 = 1 will always work
%
% Parameters that are known to work and place the solution on periodic
% orbits are:
%   v0 = 1
%   kappa = 1
%   dx = .03
%   mu = -1.1

function [u v] = full_chemIVPShoot_symmetric(v0, kappa, dx, mu)
    % Start at V_x = 0
    v(1) = v0;
    v(2) = v0+.0000391; %Determined emperically from initial condition v0 = 1
    
    % In the continuum limit u = exp(v+mu) so this should be a good guess
    % for the IC of the discrete system
    u = exp(v + mu);
    
    i=1;
    % Start the first couple calculations assuming Neumann BC at the left
    % boundary
    v(i+2) = 4*dx^2/kappa*(v(i)-u(i)) + 2*v(i) - v(i+1); % here v(i-2) = v(i+1)
    u(i+2) = 2*u(i) - u(i+1) + u(i+1)*(v(i+2)-v(i)) + u(i)*(v(i)-v(i+1));% here u(i-2) = u(i+1) and u(i-1) = u(i)
    
    i = i+1;
    
    v(i+2) = 4*dx^2/kappa*(v(i)-u(i)) + 2*v(i) - v(i-1); % here v(i-2) = v(i-1)
    u(i+2) = 2*u(i) - u(i-1) + u(i+1)*(v(i+2)-v(i)) - u(i-1)*(v(i)-v(i-1)); % here u(i-2) = u(i-1)
    
    i = i+1;
    % Stop the iteration after a maximum is reached
    while v(i+1) >= v(i-1)
%     while i<200
        % results from 0 = V_xx - V + U
        v(i+2) = 2*v(i) - v(i-2) + ...
                 4*dx^2*(v(i)-u(i))/kappa;
        
        
        % results from 0 = (U_x - U*V_x)_x
        % here we have used symmetric differenced for all of the
        % derivatives, as a result we are using a next nearest neighbor
        % approximation for the laplacian
        u(i+2) = 2*u(i) - u(i-2) + ...
                 u(i+1)*(v(i+2)-v(i)) - ...
                 u(i-1)*(v(i)-v(i-2));
        
        i = i+1;
    end
        
    %Drop the first elements of u and v
    v = v';
    u = u';
end