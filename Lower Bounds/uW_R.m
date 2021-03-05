function d = uW_R(pos,w1,w2,p)
% the function calculate Wasserstein distance on real line with respect to
% Lambda_infinity metric. pos is the union of supports of w1 and w2

%  pos, w1, w2 are 1*n arrays
%  p between 1 to infty, now only for p=1
    n = length(pos);
    pos1 = [pos(2:n) 0];
    diff = w1 - w2;
    posdif = pos1.^p - pos.^p;
    posdif(n) = [];
    %posdif(1) = [];
    
    cw = cumsum(diff);
    cw(n) = [];
    %cw(1) = [];
    
    
    d = (0.5*(abs(cw) * abs(posdif)'+ abs(diff) * (pos.^p)'))^(1/p);
end