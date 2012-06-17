function s = speed2(kappa, chi, guess)
    global L_guess
    L_guess = [0; 1];
    s = fzero(@(s) realPartLambda(s, kappa, chi), guess);
end

function r_L = realPartLambda(s, K, X)
    %find the lambda for which the wedgeProduct of the eigenvectors is zero
    global L_guess
    [L_guess,~,~] = fminsearch(@(L) norm(wedgeEigenvectors(L,s,K,X)),L_guess);
    r_L=L_guess(1);
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
    
    e_large = e_large*sign(e_large(1));
    e_small = e_small*sign(e_small(1));
    
    res = wedge(e_large,e_small);
%     sprintf('L: %f.4 + %f.4i\nwedge: %f.4 + %f.4i',real(L),imag(L),real(res),imag(res))
    res = [real(res); imag(res)]+.00001*[1; 1];
end

function res = wedge(v1,v2)
    v2(2) = -v2(2);
    v2(5) = -v2(5);
    res = sum(v1.*v2(end:-1:1));
end