%% doubleRootProblem.m
%
% Returns a function handle 'drp' (double root problem) that gives the real
% and imaginary part of the dispersion relation and it's derivative as a
% function of nu (real and imaginary), lambda (real and imaginary), speed,
% and mass.
% 
% The double roots occur at the values that make 'drp' return a vector of
% zeros.

function drp = doubleRootProblem(kappa)
    % Define the dispersion relation
    dr = @(nu,lambda,s,mass) ...
            (nu.^2+s.*nu-lambda).*...
            (kappa.*nu.^2+s.*nu-1-lambda) + ...
            mass.*nu.^2;
    
    % Define the partial derivative (with respect to nu) of the dispersion
    % relation
    dNu_dr = @(nu,lambda,s,mass) ...
            (2*nu+s).* (kappa.*nu.^2+s.*nu-1-lambda) + ...
            (nu.^2+s.*nu-lambda).*(2*kappa.*nu+s) + ...
            2*mass.*nu;
    
    % Build a problem that fsolve can handle
    drp = @(nu_r,nu_i,lambda_r,lambda_i,s,mass) ...
            [real(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             imag(dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             real(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass));
             imag(dNu_dr(nu_r+1i*nu_i,lambda_r+1i*lambda_i,s,mass))];

end