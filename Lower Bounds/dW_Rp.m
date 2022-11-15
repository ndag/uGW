function cost=dW_Rp(pos1, pos2, mu1, mu2, p)

%   A Matlab function that implements the Wasserstein distance of order p on the real line.

%   pos1  - support of mu1
%   pos2  - support of mu2
%   mu1 - probability vector
%   mu2 - probability vector
%   p   - real number >=1

%   
% Returns:
%   res   - the Wasserstein distance between the probability measures mu1 and mu2 to the power of p.
    [pos1,tmp1] =sort(pos1);
    mu1=mu1(tmp1);

    [pos2,tmp2] =sort(pos2);
    mu2=mu2(tmp2);
    
    cost = 0;
    
    n = length(mu1);
    m = length(mu2);

    i = 1;
    w_i = mu1(1);
    j = 1;
    w_j = mu2(1);

    m_ij = 0;
    
    while(true)
        m_ij = abs(pos1(i) - pos2(j))^p;
        if (w_i < w_j || j == m)
            cost = cost+m_ij * w_i;
            i = i+1;
            if (i == n+1)
                break;
            end    
            w_j = w_j-w_i;
            w_i = mu1(i);
        else
            cost = cost+ m_ij * w_j;
            j = j+1;
            if (j == m+1)
                break;
            end
            w_i =w_i- w_j;
            w_j = mu2(j);
        end    
    end
   end