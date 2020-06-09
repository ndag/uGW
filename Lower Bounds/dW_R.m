function d = dW_R(pos,w1,w2)
% the function calculate Wasserstein distance on real line with respect to
% Delta_1 metric. pos is the union of supports of w1 and w2

%  pos, w1, w2 are 1*n arrays
%  p between 1 to infty, now only for p=1
    n = length(pos);
    pos1 = [pos(2:n) 0];
    diff = w1 - w2;
    posdif = pos1 - pos;
    posdif(n) = [];
    %posdif(1) = [];
    
    cw = cumsum(diff);
    cw(n) = [];
    %cw(1) = [];
    
    
    d = abs(cw) * abs(posdif)';
end