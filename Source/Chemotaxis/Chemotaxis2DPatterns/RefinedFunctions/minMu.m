function mu = minMu(kappa,b,mass,T)
    temp = (kappa.*(2*pi./T).^2+1)./b;
    mu = min(log(mass)-b.*mass,log(temp)-b.*temp);
end