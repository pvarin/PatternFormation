%% chemBVP.m
%
% Defines the boundary value problem as a function in N points using finite
% differences
%   0 = u_xx - (u*v_x)_x
%   0 = kappa*v_xx - v + u
%
% We are considering periodic solutions on length L with N points.
% We are using the symmetric derivative approximation and 

function bvp = chemBVP(N)
    %% Construct the Derivative and Laplacian Matrices
    Lap = laplacian(1,N,'periodic');
    D = derivative(1,N,'symmetric','periodic');

    %% Compute the Error
    u = @(sol) sol(1:N);
    v = @(sol)sol(N+1:2*N);
    dx = @(L) L/N;
    
    bvp = @(sol,L,kappa,mass) [(Lap/dx(L)^2*u(sol)-D*(u(sol).*(D*v(sol)))/dx(L)^2)+sol(end-1);...
                          (kappa*Lap/dx(L)^2*v(sol)+u(sol)-v(sol));
                          mean(u(sol))-mass;
                          sol(1)-sol(N)] + ...
                          sol(end)*([sol(1:2*N);0;0]-[sol(N);sol(1:N-1);sol(2*N);sol(N+1:2*N-1);0;0]);
end