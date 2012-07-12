function mainTest
    %% Constants and Parameters
    % Fundamental Constants
    global e
    e=exp(1);
    
    % Model Parameters
    kappa = 1;
    b = 1;
    mass = 1.01; % unused as of now
    LDes = 20; % unused as of now
    
    % Secondary Parameters
    maxMu = maximumMu(kappa,b,mass,LDes);
    mu = maxMu-.1;
    v0 = 1/b;
    dx = 0.01;
    
    %% Create an Initial Guess
    % Find a solution to the differential equation
    [v N L] = initialShoot(v0,dx);
    plot(v);drawnow;shg
    bvp = chemBvp(kappa,b,N);
    v = chemSolve(v,bvp);
    plot(v);drawnow
    
    %% Body
    % Follow the solution in L in increments of ~0.1
    fprintf('Desired Period: %f\n',LDes)
    Lstep = (LDes-L)/floor((LDes-L)/1);
    while abs(LDes-L) > 1e-3
        fprintf('    Current Period: %f\n',L)
        L = L+Lstep;
        v = chemSolve(v,bvp);
        plot(v);drawnow
    end
    fprintf('Final Period: %f\n',L)
    
    % Follow the solution in mu until the correct mass is acquired
    % Note: increasing mu will decrease the mass
    fprintf('Desired Mass: %f \n',mass)
    muStep = .1*sign(calcMass()-mass);
    
    while abs(mass-calcMass()) > 1e-5
        fprintf('    Current Mass: %f \n',calcMass())
        if (mu + muStep > maxMu)
            warning('Mu is greater than the max!!!')
        end
        mu = mu + muStep;
        v = chemSolve(v,bvp);
        plot(v);drawnow
        
        if sign(calcMass()-mass)*muStep < 0
            mu = mu - muStep;
            muStep = .5*muStep;
        end
    end
    fprintf('Final Mass: %f \n',calcMass())
    plot(v)
    
    % Refine the grid
    v = interp(v,5);
    N = length(v);
    dx = L/N;

    bvp = chemBvp(kappa,b,N);
    v=chemSolve(v,bvp);
    plot(v);drawnow
    
    %% Nested Function Definitions
    function [v N L]=initialShoot(v0,dx)
        %{
        Solves the IVP 
        Used to generate an inital guess for V
        %}
        v(1)=v0;
        v(2)=v(1);
        j=2;
        while (v(j) >=v(j-1))
            v(j+1)=2*v(j)-v(j-1)+dx^2/kappa*(-exp(mu+v(j))+v(j));
            j=j+1;
        end
        N=j-1;
        L=N*dx*2;
        v=v(2:length(v))';
    end

    function vsol=chemSolve(initial,bvp)
        %{
        Finds the solutions to the boundary value problem defined in bvp
        using fsolve
        %}

        options = optimset('TolX',1e-8,'Display','off');
        vsol=fsolve(@(v) bvp(v,mu,L),initial,options);
    end

    function m=calcMass()
        m = mean(v);
    end
end

%% Helper Function Definitions
function bvp=chemBvp(kappa,b,N)
        %{
        Defines boundary value problem as function in N points using finite differences
        kappa v_xx+ exp(mu+v)-v=0

        We are solving for even periodic solutions
        This function is defined by N points on its half period (L/2)
        This is equivalent to the Neumann problem on the full period (L)
        %}

        E=ones(N,1);
        Lap=spdiags([E -2*E E],-1:1,N,N);
        Lap(1,1)=-1; Lap(N,N)=-1;

        % k*V_xx = V - exp(u+b*V)
        % dx = L/(2*N)
        %}
        bvp = @(v,mu,L) kappa/(L/(2*N))^2*Lap*v + exp(mu+b*v) - v;
end

function mu = maximumMu(kappa,b,mass,T)
    temp = (kappa.*(2*pi./T).^2+1)./b;
    mu = min(log(mass)-b.*mass,log(temp)-b.*temp);
end