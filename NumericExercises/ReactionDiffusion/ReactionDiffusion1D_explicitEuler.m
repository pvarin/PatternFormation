% model parameters
a = .5;
kappa = .1;
gamma = 1;

% similation parameters
dt = .00001;
dx = .1;
numsteps = 10000;
N = 100;

% initial conditions
C0 = zeros(1,N);
E0 = a*ones(1,N); E0(1,1) = E0(1,1) + .01;

% initialize the matrices
C = zeros(numsteps,N); E = zeros(numsteps,N);
C(1,:) = C0; E(1,:) = E0;

%Explicit Euler Method
for i=2:numsteps
    % domain interior
    n = 2:N-1;
    f = E(i-1,n).*(1-E(i-1,n)).*(E(i-1,n)-a)+gamma*C(i-1,n);
    C(i,n) = C(i-1,n)-(C(i-1,n-1)-2*C(i-1,n)+C(i-1,n+1))*dt/dx^2-f*dt;
    E(i,n) = E(i-1,n)-(E(i-1,n-1)-2*E(i-1,n)+E(i-1,n+1))*dt/dx^2+f*dt;
    
    % domain boundaries (Neumann Conditions)
    n = 1;
    f = E(i-1,n).*(1-E(i-1,n)).*(E(i-1,n)-a)+gamma*C(i-1,n);
    C(i,n) = C(i-1,n)-(-C(i-1,n)+C(i-1,n+1))*dt/dx^2-f*dt;
    E(i,n) = E(i-1,n)-(-E(i-1,n)+E(i-1,n+1))*dt/dx^2+f*dt;
    
    n = N;
    f = E(i-1,n).*(1-E(i-1,n)).*(E(i-1,n)-a)+gamma*C(i-1,n);
    C(i,n) = C(i-1,n)-(C(i-1,n-1)-C(i-1,n))*dt/dx^2-f*dt;
    E(i,n) = E(i-1,n)-(E(i-1,n-1)-E(i-1,n))*dt/dx^2+f*dt;
end

figure(1)
pcolor(linspace(0,dx*N,N),linspace(0,dt*numsteps,numsteps),C)
shading flat
drawnow

figure(2)
hold all
for i=1:10:numsteps
    clf
    shg
    plot(linspace(0,dx*N,N),C(i,:))
    hold all
    plot(linspace(0,dx*N,N),E(i,:))
    axis([0 dx*N -2 2])
    drawnow
end