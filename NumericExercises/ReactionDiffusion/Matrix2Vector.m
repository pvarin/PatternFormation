function V = Matrix2Vector(M)
    % extract matrix size
    s = size(M);
    n = s(1);
    m = s(2);
    
    % initialize vector
    V = zeros(n*m,1);
    
    for i=1:m
        V(((i-1)*n+1):n*i) = M(:,i);
    end
end