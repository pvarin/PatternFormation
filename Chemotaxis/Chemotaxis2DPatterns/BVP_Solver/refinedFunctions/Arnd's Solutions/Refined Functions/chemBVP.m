%% chemBVP.m
%
% Defines the boundary value problem as a function in N points using finite
% differences
%   0 = u_xx - (u*v_x)_x
%   0 = kappa*v_xx - v + u
%
% We are considering periodic solutions on length L with N points.
% We are using the symmetric derivative approximation and 

function sol = chemBVP(sol,L,kappa)
    %% Construct the Derivative and Laplacian Matrices
    N = length(sol)/2;
    dx=L/N;
    Lap = laplacian(dx,N,'periodic');
    D = derivative(dx,N,'symmetric','periodic');

    %% Compute the Error
    u=sol(1:N);
    v=sol(N+1:2*N);

    sol=[(Lap*u-D*(u.*(D*v)));...
         (kappa*Lap*v+u-v)];
end