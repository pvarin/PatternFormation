function [P_,D_]=sortem(P,D)

D_=diag(sort(diag(D),'descend')); % make diagonal matrix out of sorted diagonal values of input D
[~, ind]=sort(diag(D),'descend'); % store the indices of which columns the sorted eigenvalues come from
P_=P(:,ind); % arrange the columns in this order

end