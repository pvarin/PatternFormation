%% full_chemIVPShoot.m
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

function [u v] = full_chemIVPShoot(v0, kappa, dx, mu)
    % Start at V_x = 0
    v(1) = v0;
    v(2) = v0;
    
    % In the continuum limit u = exp(v+mu) so this should be a good guess
    % for the IC of the discrete system
    u = exp(v + mu);
    
    % Stop the iteration after a maximum is reached
    i=2;
    while v(i) >= v(i-1)
        % results from 0 = V_xx - V + U
        v(i+1) = 2*v(i) - v(i-1) + ...
                 dx^2*(v(i)-u(i))/kappa;
        
        
        % results from 0 = (U_x - U*V_x)_x
        % here we have used the right difference for the inner derivatives
        % and the left difference for the outer derivatives
        %   Note that this is also consistent with the approximation for
        %   the laplacian in the V equation
        u(i+1) = 2*u(i) - u(i-1) + ...
                 (u(i)-u(i-1))*(v(i)-v(i-1)) + ...
                 u(i)*(v(i-1)-2*v(i)+v(i+1));
        
        i = i+1;
    end
        
    %Drop the first elements of u and v
    v = v(2:end)';
    u = u(2:end)';
end