function res = wedge(v1,v2)
    v2(2) = -v2(2);
    v2(5) = -v2(5);
    res = sum(v1.*v2(end:-1:1));
end