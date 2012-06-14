function V = fixedPoint2(mu,b)
    equilCond = @(V_) V_-exp(mu+b*V_);
    V = fzero(equilCond,[1/b,-2*(mu-1/b)/b]);
end