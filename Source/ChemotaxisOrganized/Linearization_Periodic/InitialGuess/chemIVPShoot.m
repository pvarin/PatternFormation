%% chemIVPShoot.m
%
% Uses shooting to develop a solution, over a full period, to the
% chemotaxis equations
%
% This function is used to build an initial guess for the boundary value
% problem for the chemotaxis equations.
% 
% Parameter values to create a periodic solution (determined emperically):
%   v0:     e+.8
%   dx:     0.0060016
%   mu:     1-e
%   kappa:  1

function [v N L] = chemIVPShoot(v0,dx,mu,kappa)
    %% Initial Conditions
    % Begin with homogeneous Neumann on the boundary
    %   i.e. v(1) = v(2)
    v(1)=v0;
    v(2)=v(1);
    
    %% Build solution
    % Iterate until the solution begins to increase
    j=2;
    while (v(j)-v(j-1) <=0)
        v(j+1)=2*v(j)-v(j-1)+dx^2/kappa*(-exp(mu+v(j))+v(j));
        j=j+1;
    end
    % Iterate again until the solution begins to decrease
    while (v(j)-v(j-1) >=0)
        v(j+1)=2*v(j)-v(j-1)+dx^2/kappa*(-exp(mu+v(j))+v(j));
        j=j+1;
    end
    
    %% Return Results
    % Drop the first and last points
    N=j-2;
    L=N*dx;
    v=v(2:end-1)';
end