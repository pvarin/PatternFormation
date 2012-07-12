function [u mu] = chemSolveMass(u_guess,mu_guess,L,mass,bvp)
    opts = optimset('TolX',1e-8,'Display','off');
    guess = [u_guess; mu_guess];
    sol = fsolve(@(sol_) bvp(sol_,mass,L),guess,opts);
    u = sol(1:end-1);
    mu = sol(end);
end