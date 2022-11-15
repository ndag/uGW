function d = dgwslb(dx, dy, mux, muy, p)

%   A Matlab function for computing the second lower bound of dGW of two metric measure spaces.

%   dx  - ultrametric distance matrix
%   dy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res - pth-uGW-second lower bound between dx and dy with probability measures mux and muy, respectively.

    [u1_dist,prob1]= distance_distribution(dx,mux); 
    [u2_dist,prob2]= distance_distribution(dy,muy);
    
   
    
    d = (dW_Rp(u1_dist,u2_dist,prob1,prob2,p))^(1/p);
    
end