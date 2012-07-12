kappas = .1:.2:10;
for i=1:length(kappas)
    kappa = kappas(i);
    mainOnePeriod
    filename = sprintf('../../../../initialSolution_kappa_%d',cast(kappa*10,'uint8'));
    save(filename,'kappa','v','u','L')
end