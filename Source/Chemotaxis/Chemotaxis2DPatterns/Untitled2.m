func = @(x_,mu_) exp(mu_+x_)-x_;
% func_approx = @(x_,mu_) 
for mu = linspace(-1.1,-10)
    x = fzero(@(x_) func(x_,mu),[1,100]);
    plot(-mu,x,'r.')
    hold on
    drawnow
end
mu = linspace(-1.1,-10);
plot(-mu,-mu+5)
% plot(mu,abs((1+sqrt(1-2*exp(2*mu+2)))./exp(mu+1)))