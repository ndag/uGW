function d = uW_R(pos, mu1, mu2, p)

%   This function implements the Wasserstein distance on (R,Delta_infinity).

%   pos - joint support of mu1 and mu2
%   mu1 - probability vector that indicates the mass of mu1 on the corresponding position in the joint support
%   mu2 - probability vector that indicates the mass of mu2 on the corresponding position in the joint support
%   p   - real number >=1

%   
% Returns:
%   res   - the (R,Delta_infinity)-Wasserstein distance between the probability measures mu1 and mu2 to the power of p.

n = length(pos);
    pos1 = [pos(2:n) 0];
    diff = mu1 - mu2;
    posdif = pos1.^p - pos.^p;
    posdif(n) = [];
    %posdif(1) = [];
    
    cw = cumsum(diff);
    cw(n) = [];
    
    
    d = (0.5*(abs(cw) * abs(posdif)'+ abs(diff) * (pos.^p)'));
end