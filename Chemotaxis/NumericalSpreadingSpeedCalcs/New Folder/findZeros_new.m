function res = findZeros_new(func,varargin)

    x = 0.01; y = 0.01;
    
    switch nargin
        case 1
            %find a starting place
            x = abs(fminsearch(@(a) -func(abs(a),0),0.01));
            y = abs(fminsearch(@(b) func(x,abs(b)),0.01));
            plot3(y,x,func(x,y),'g*')
            condition = true
            c = fminsearch(@(c) func(x*c*1.01,y*c),1);
            if c>0
                x = 1.01*c*x;
                y = c*y;
            end
            c = fminsearch(@(c) func(x*c,y*c*1.01),1);
            if c>0
                x = c*x;
                y = 1.01*c*y;
            end
            vect1 = [1 0];
            vect2 = [0 1];
            plot3(y,x,func(x,y),'r*');
            iter = 0;
        case 2
            %extract the guess
            x = varargin{1}(1);
            y = varargin{1}(2);
            vect1 = [1 0];
            vect2 = [0 1];
        case 3
            x = varargin{1}(1);
            y = varargin{1}(2);
            vect1 = [1 0];
            vect2 = [0 1];
            warning('findZeros: third argument ignored')
        case 5
            %extract the guess and the directions
            x = varargin{1}(1);
            y = varargin{1}(2);
            vect1 = varargin{2};
            vect2 = varargin{3};
            iter = varargin{4};
        otherwise
            warning('findZeros: extra arguments ignored')
    end
    
    data1 = zeros(2);
    data2 = zeros(2);
    for i=1:2
        x = fminsearch(@(a) func(vect1(1)*a+x,vect1(2)*a+y),0)*vect1(1)+x;
        data1(i,:) = [x y];
        plot3(y,x,func(x,y),'r*'); drawnow
        
        y = fminsearch(@(b) func(vect2(1)*b+x,vect2(2)*b+y),0)*vect2(2)+y;
        data2(i,:) = [x y];
        plot3(y,x,func(x,y),'r*'); drawnow
    end
    plot3(y,x,func(x,y),'g*'); drawnow
    
    if (func(x,y)<1e-5)
        res = y;
    else
        vect1 = data1(2,:)-data1(1,:);
        vect2 = data2(2,:)-data2(1,:);
        res = findZeros_new(func,[x,y],vect2,vect1,iter+1);
    end
end