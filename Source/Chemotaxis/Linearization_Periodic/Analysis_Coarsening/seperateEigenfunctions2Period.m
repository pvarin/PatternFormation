load('../../Data/equilibriumSolutions_kappa_10_mass_11_50.mat')

U = [U; U(end:-1:1,:)];
V = [V; V(end:-1:1,:)];
periods = periods*2;

N = length(U(:,1));
numSolutions = length(U(1,:));
numEigen = 4;

last_eVects = zeros(2*N,numEigen,2);
all_eVals = zeros(numSolutions,numEigen);
all_eVects = zeros(2*N,numEigen,numSolutions);

for i=1:numSolutions
    Lin = chemLinear(U(:,i),V(:,i),periods(i),kappa);
    [eVects eVals] = eigs(Lin,numEigen,3);
    [eVects eVals] = sortEigen(eVects, eVals);
    
    if i==1
        last_eVects(:,:,1) = eVects;
        last_eVects(:,:,2) = eVects;
        all_eVals(i,:) = eVals;
        all_eVects(:,:,i) = eVects;
        continue
    end
    
%     plot(eVects)
%     
    index = zeros(numEigen,1);
    for j=1:numEigen
        correlation = 0;
        
        for k=1:numEigen
            if correlation < abs(last_eVects(:,j,1)'*eVects(:,k)) + ...
               abs(last_eVects(:,j,2)'*eVects(:,k))
                
                correlation = abs(last_eVects(:,j,1)'*eVects(:,k)) + ...
                abs(last_eVects(:,j,2)'*eVects(:,k));
                index(j) = k;
            end
        end
    end
    
    % update last_eVects
    last_eVects(:,:,1:end-1) = last_eVects(:,:,2:end);
    last_eVects(:,:,end) = eVects(:,index);
    
    % log data
    all_eVals(i,:) = eVals(index);
    all_eVects(:,:,i) = eVects(:,index);
    i
end

plot(masses,all_eVals);
save('Instability2Periods.mat','all_eVects','all_eVals','periods','U','V','masses','kappa')