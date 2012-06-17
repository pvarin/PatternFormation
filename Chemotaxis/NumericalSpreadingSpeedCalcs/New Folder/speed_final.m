function [T s] = findSpeedPeriod(kappa, chi)
    % Minimize the modulus of the wedge product for a given kappa and chi
    func = @(a,b) eigenvectors(a*1i,b,kappa,chi);
    r = findZeros(func);
    
    % Extract the speed and the selected period
    s = r(2);
    T = 2*pi*s/r(1);
end

function res = eigenvectors(L, s, K, X)
    %Construct the Matrix
    A2 = [[-s       X*(L+1)/K   -X*s/K      0       0           0];
          [0        0           1           1       0           0];
          [0        (L+1)/K     -s/K        0       1           0];
          [0        L-X/K       0           -s      1           X*s/K];
          [1/K      0           L-X/K       (L+1)/K -(s+s/K)    X*(L+1)/K];
          [0        1/K         0           0       0           -s/K]];
    % Compute the eigenvectors of the largest and smallest eigenvalues
    % (by real part)
    [e_large, ~] = eigs(A2,1,'lr');
    [e_small, ~] = eigs(A2,1,'sr');
    
    % Return the modulus of the wedge product of these eigenvectors
    res = abs(wedge(e_large,e_small));
end

function res = wedge(v1,v2)
    % Compute the wedge product of v1 and v2
    v2(2) = -v2(2);
    v2(5) = -v2(5);
    res = sum(v1.*v2(end:-1:1));
end

function res = findZeros(func,varargin)

    x = 0.01; y = 0.01;
    
    switch nargin
        case 1 % This is the first time the function is called
            % Find a starting place
            x = abs(fminsearch(@(a) -func(abs(a),0),0.01));
            y = abs(fminsearch(@(b) func(x,abs(b)),0.01));
            
            % Caveat 1
            c = fminsearch(@(c) func(x*c*1.01,y*c),1);
            if c>0
                x = 1.01*c*x;
                y = c*y;
            end
            % Caveat 2
            c = fminsearch(@(c) func(x*c,y*c*1.01),1);
            if c>0
                x = c*x;
                y = 1.01*c*y;
            end
            
            % Define initial directions in which to minimize
            vect1 = [1 0];
            vect2 = [0 1];
        case 3
            % Extract the guess
            x = varargin{1}(1);
            y = varargin{1}(2);
            % Extract the directions
            vect1 = varargin{2};
            vect2 = [vect1(2), -vect1(1)];% Orthogonal to vect1
    end
    
    data1 = zeros(2);
    data2 = zeros(2);
    for i=1:2
        % Minimize in the direction of vect1
        x = fminsearch(@(a) func(vect1(1)*a+x,vect1(2)*a+y),0)*vect1(1)+x;
        data1(i,:) = [x y];% store data
        
        % Minimize in the direction of vect2
        y = fminsearch(@(b) func(vect2(1)*b+x,vect2(2)*b+y),0)*vect2(2)+y;
        data2(i,:) = [x y];% store data
    end
    
    if (func(x,y)<1e-5)
        % If the result is within the tolerance stop recursion
        res = [x, y];
    else
        % Determine new directions to minimize in given the data
        vect1 = data1(2,:)-data1(1,:);
        vect2 = data2(2,:)-data2(1,:);
        res = findZeros(func,[x,y],(vect1+vect2)/2);
    end
end