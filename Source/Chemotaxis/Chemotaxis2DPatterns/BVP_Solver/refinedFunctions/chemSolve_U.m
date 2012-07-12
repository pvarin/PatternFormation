function u = chemSolve_U(u,mu,L,bvp)
    opts = optimset('TolX',1e-8,'Display','off');
    u = fsolve(@(u) bvp(u,mu,L),u,opts);
end