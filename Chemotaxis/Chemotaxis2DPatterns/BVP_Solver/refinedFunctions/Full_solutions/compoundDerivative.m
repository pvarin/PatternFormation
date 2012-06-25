function D = compoundDerivative(dx,N,style,numPeriods)
    % Periodic boundary conditions
    % N is the number of points in the half period
    
    D_temp1 = derivative(dx,N,style,'');
    D_temp2 = -D_temp1(end:-1:1,end:-1:1);
    D = zeros(N*2*numPeriods,N*2*numPeriods);
    
    for k=1:numPeriods
        i = (k-1)*N*2+1:(2*k-1)*N;
        
        D(i,i) = D_temp1;
        D(i+N,i+N) = D_temp2;
        
        if strcmp(upper(style),'RIGHT')
            D(i(end),i(end)+1) = 1/dx;
            D(i(end)+1,i(end)) = -1/dx;
        end
        if strcmp(upper(style),'SYMMETRIC')
            D(i(end),i(end)+1) = 1/(dx*2);
            D(i(end)+1,i(end)) = -1/(dx*2);
        end
    end
    

    if strcmp(upper(style),'LEFT')
        D(1,end) = -1/dx;
        D(end,1) = 1/dx;
    end
    if strcmp(upper(style),'SYMMETRIC')
        D(1,end) = -1/(dx*2);
        D(end,1) = 1/(dx*2);
    end
end