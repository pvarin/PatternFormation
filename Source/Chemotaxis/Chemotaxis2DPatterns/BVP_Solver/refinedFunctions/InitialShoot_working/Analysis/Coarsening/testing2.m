load('Instability1Period.mat')

% find the stable eigenvectors for period 1
u = U(:,100);
v = V(:,100);
L = periods(100);
N =length(u);

stableU = all_eVects(1:end/2,1,100);
stableV = all_eVects(end/2+1:end,1,100);
% stableU_ = u([2:N,1])-u([N,1:N-1]);
% stableV_ = v([2:N,1])-v([N,1:N-1]);
% n = norm([stableU_; stableV_]);
% stableU_ = stableU_/n;
% stableV_ = stableV_/n;
% 
% plot(stableU+stableU_)
% hold all
% plot(stableV+stableV_)
% return

Lin = chemLinear(u,v,L,kappa);

% sol = fsolve(@(sol) Lin*sol,[stableU; stableV]);
% stableU = sol(1:N);
% stableV = sol(N+1:2*N);

figure(2)
plot([stableU stableV])

% find the stable and unstable eigenvectors for period 2
eigU2 = [stableU; stableU(end:-1:1)];
eigV2 = [stableV; stableV(end:-1:1)];

u2 = [u; u(end:-1:1)];
v2 = [v; v(end:-1:1)];
L = L*2;
N = N*2;
dx = L/N;

x = dx/2:dx:L;
x = x';
gamma = 2*pi/L;
phase = [exp(1i*gamma*x); exp(1i*gamma*x)];
stableEig2 = [eigU2; eigV2];
% stableEig2 = fsolve(@(sol) Lin2*sol, stableEig2);
unstableEig2 = stableEig2.*phase;
% unstableEig2 = unstableEig2;

% check that these are eigen vectors of the 2-period system
Lin2 = chemLinear(u2,v2,L,1);
figure(3)
plot(Lin2*stableEig2)
figure(4)
plot(real([(Lin2*unstableEig2) unstableEig2]))
figure(5)
plot((Lin2*unstableEig2)./unstableEig2)