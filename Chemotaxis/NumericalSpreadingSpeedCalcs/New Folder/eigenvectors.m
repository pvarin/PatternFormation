function res = eigenvectors(L, s, K, X)
    A2 = [[-s       X*(L+1)/K   -X*s/K      0       0           0];
          [0        0           1           1       0           0];
          [0        (L+1)/K     -s/K        0       1           0];
          [0        L-X/K       0           -s      1           X*s/K];
          [1/K      0           L-X/K       (L+1)/K -(s+s/K)    X*(L+1)/K];
          [0        1/K         0           0       0           -s/K]];
    [e_large, ev_large] = eigs(A2,1,'lr');%eigenvector corresponding to the eigenvector with the largest real part
    [e_small, ev_small] = eigs(A2,1,'sr');%eigenvector corresponding to the eigenvector with the smallest real part
    
    e_large = e_large*sign(e_large(1));
    e_small = e_small*sign(e_small(1));
    
    res = abs(wedge(e_large,e_small));
%     res = sign(real(res))*abs(res);
end