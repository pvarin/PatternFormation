function rhs=chembvp(v,N,mu,L,kappa)

%defines boundary value problem as function in N points using finite differences
% kappa v_xx+ exp(mu+v)-v=0

% note that we are thinking of a periodic problem on length 2L with N points, even functions, which is equivalent to a Neumann problem on length L domain

dx=L/N/2;
e=ones(N,1);
Lap=spdiags([e -2*e e],-1:1,N,N);
Lap(1,1)=-1;Lap(N,N)=-1;
Lap=kappa/dx^2*Lap;

%  rhs=Lap*v+exp(mu+v)-v;%  
rhs=Lap*v+exp(mu+v)-v;

end


