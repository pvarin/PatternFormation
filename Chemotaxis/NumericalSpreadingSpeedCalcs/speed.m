function s = speed(kappa, chi, guess)
    global L_guess
    L_guess = [1; 1];
    s = fzero(@(s) realPartLambda(s, kappa, chi), guess);
end

function r_L = realPartLambda(s, K, X)
    %find the lambda for which the wedgeProduct of the eigenvectors is zero
    global L_guess
    L_guess = fsolve(@(L) wedgeEigenvectors(L,s,K,X),L_guess);
    r_L=L(1);
end

function res = wedgeEigenvectors(L, s, K, X)
    L = L(1) + 1i*L(2);
    A2 = [[-s       X*(L+1)/K   -X*s/K      0       0           0];
          [0        0           1           1       0           0];
          [0        (L+1)/K     -s/K        0       1           0];
          [0        L-X/K       0           -s      1           X*s/K];
          [1/K      0           L-X/K       (L+1)/K -(s+s/K)    X*(L+1)/K];
          [0        1/K         0           0       0           -s/K]];
    [e_large, ~] = eigs(A2,1,'lr');%eigenvector corresponding to the eigenvector with the largest real part
    [e_small, ~] = eigs(A2,1,'sr');%eigenvector corresponding to the eigenvector with the smallest real part
    
    res = wedge(e_large,e_small);
    res = [real(res); imag(res)];%FIXME
end

function res = wedge(v1,v2)
    v2(2) = -v2(2);
    v2(5) = -v2(5);
    res = sum(v1.*v2(end:-1:1));
end