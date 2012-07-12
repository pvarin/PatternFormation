kappa = 1;
b = 1;
mu = linspace(-1-log(b)-.01,-10,20);
mu_ = mu;
V_fixed = @(V_,i_) exp(mu_(i_)+b*V_)-V_;

Vmin = fzero(@(V_) V_fixed(V_,length(mu)),[0,1/b]);
Vmax = fzero(@(V_) V_fixed(V_,length(mu)),[1/b,100]);%FIXME: this is a bad upper bound

V0 = linspace(Vmin,Vmax,20);

[mu, V0] = meshgrid(mu,V0);

period = zeros(size(mu));
mass = zeros(size(mu));

for i=1:length(mu(1,:))
    i
    Vmin = fzero(@(V_) V_fixed(V_,i),[0,1/b]);
    Vmax = fzero(@(V_) V_fixed(V_,i),[1/b,100]);%FIXME: this is a bad upper bound 
    for j=1:length(V0(:,1))
        if (V0(i,j) <= Vmin) || (V0(i,j) >= Vmax)
            continue
        end
        [period(i,j) mass(i,j)] = IVP_onePeriod(kappa, b, mu(i,j), V0(i,j));
    end
    
end
surf(mu(1,:),V0(:,1),period);
surf(mu(1,:),V0(:,1),mass);