load('equilibriumSolutions_kappa_10_mass_11_50.mat')

index = [130];
numsteps = [120000];
filenames = {'CoarseningChemotaxisData_mass_2-8788.mat'};         
for i=1:length(index)
    % create initial conditions
    U0 = [U(:,index(i)); U(end:-1:1,index(i))];
    V0 = [V(:,index(i)); V(end:-1:1,index(i))];
    N = length(U0);
    dx = 2*periods(index(i))/N;
    
    %perturb the initial conditions
    dU0_dm = [U(:,index(i)-1); U(end:-1:1,index(i)-1)];
    dV0_dm = [V(:,index(i)-1); V(end:-1:1,index(i)-1)];
    dU0_dm = U0 - dU0_dm; dU0_dm = dU0_dm/norm(dU0_dm);
    dV0_dm = V0 - dV0_dm; dV0_dm = dV0_dm/norm(dV0_dm);
    
    dU0_dx = U0([2:N 1])-U0([N 1:N-1]);
    dV0_dx = V0([2:N 1])-V0([N 1:N-1]);
    dU0_dx = dU0_dx/norm(dU0_dx);
    dV0_dx = dV0_dx/norm(dV0_dx);
    
    epsU = (dU0_dm + dU0_dx).*cos(linspace(0,2*pi,N))';
    epsV = (dV0_dm + dV0_dx).*cos(linspace(0,2*pi,N))';
    
    U0 = U0 + epsU;
    V0 = V0 + epsV;
    
    %shift the solutions so that the peaks are not on the boundary    %shift the IC so that the peaks will coalesce in the center
    U0 = U0([N/4+1:N,1:N/4]);
    V0 = V0([N/4+1:N,1:N/4]);
    
    Chemotaxis1D('u0',U0,'v0',V0,'dt',.01,...
                 'numsteps',numsteps(i),'dx',dx,...
                 'savedata',filenames{i})
end