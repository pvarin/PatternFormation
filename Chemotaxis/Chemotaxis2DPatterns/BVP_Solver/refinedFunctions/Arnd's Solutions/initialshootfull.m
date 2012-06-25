function [v N L]=initialshootfull(v0,dx,mu,kappa)
    v(1)=v0;
    v(2)=v(1);
    j=2;
    while (v(j)-v(j-1) <=0)
        v(j+1)=2*v(j)-v(j-1)+dx^2/kappa*(-exp(mu+v(j))+v(j));
        j=j+1;
    end
    while (v(j)-v(j-1) >=0)
        v(j+1)=2*v(j)-v(j-1)+dx^2/kappa*(-exp(mu+v(j))+v(j));
        j=j+1;
    end
    N=j-2;
    L=N*dx;
    v=v(2:end-1)';
end