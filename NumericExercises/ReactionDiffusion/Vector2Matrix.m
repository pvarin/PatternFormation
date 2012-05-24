function M = Vector2Matrix(V,m,n)
    % initialize the matrix
    M = zeros(n,m);
    %construct the matrix
    for i=1:m
        M(:,i) = V(((i-1)*n+1):i*n);
    end
end