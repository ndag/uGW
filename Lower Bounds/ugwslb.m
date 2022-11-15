function d = ugwslb(ux, uy, mux, muy, p)

%   A function for computing the second lower bound of uGW of two ultrametric measure spaces.

%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res - pth-uGW-second lower bound between ux and uy with probability measures mux and muy, respectively.
    
    [u1_dist,prob1]= distance_distribution(ux,mux); 
    [u2_dist,prob2]= distance_distribution(uy,muy);
    
    %generate position by taking union
    pos = union(u1_dist,u2_dist);
    pos = pos';
    

    n = length(pos);
    w1 = zeros(1,n);
    w2 = zeros(1,n);
    
    j = 1;
    k = 1;
    for i = 1:n
       if j<=length(u1_dist) && pos(i)==u1_dist(j)
           w1(i)=prob1(j);
           j=j+1;
       end
      if k<=length(u2_dist) && pos(i)==u2_dist(k)
           w2(i)=prob2(k);
           k=k+1;
      end
    end
    
    d = uW_R(pos,w1,w2,p)^(1/p);
    
end