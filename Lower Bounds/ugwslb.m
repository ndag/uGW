function d = ugwslb(u1,u2,p)
%the function calculates the second lower bound of ugw between two
%ultrametric spaces with uniform measure
    n1 = length(u1);
    n2 = length(u2);
    v1 = reshape(u1,[],1);
    v2 = reshape(u2,[],1);
    [GC1,GR1] = groupcounts(v1);
    [GC2,GR2] = groupcounts(v2);

    %generate position by taking union
    pos = union(GR1,GR2);
    pos = pos';
    
    %generate weight; now use uniform measure

    n = length(pos);
    w1 = zeros(1,n);
    w2 = zeros(1,n);

    for i = 1: length(GC1)
        w1(pos == GR1(i)) = GC1(i);
    end

    for i = 1: length(GC2)
        w2(pos == GR2(i)) = GC2(i);
    end
    
    w1 = w1 ./ (n1 * n1);
    w2 = w2 ./ (n2 * n2);
    
    d = uW_R(pos,w1,w2,p);
    
end