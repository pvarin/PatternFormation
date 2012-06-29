function [eVect eVals] = sortEigen(eVect,eVals)
    [eVals index] = sort(diag(eVals),'descend');
    eVect = eVect(:,index);
end