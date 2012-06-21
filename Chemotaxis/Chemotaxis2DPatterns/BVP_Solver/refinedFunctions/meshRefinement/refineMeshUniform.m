function V = refineMeshUniform(kappa,V,L,meshFactor)
    mass = mean(V);
    mu = log(mass)-log(mean(exp(V)));
    V = interp1(V,meshFactor);
    N = length(V);
    
    bvp = chemBvpMass(kappa,N);
    
    function [v mu] = chemSolveMass(v_guess,mu_guess,L,mass,bvp)
        opts = optimset('TolX',1e-8,'Display','off');
        guess = [v_guess; mu_guess];
        sol = fsolve(@(sol_) bvp(sol_,mass,L),guess,opts);
        v = sol(1:end-1);
        mu = sol(end);
    end

    [V ~] = chemSolveMass(V, mu, L, mass, bvp);
end