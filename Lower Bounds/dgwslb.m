function d = dgwslb(u1,u2,mu1,mu2, p)
% the function calculates the second lower bound of dgw between two
% ultrametric spaces with arbitrary measures

    [u1_dist,prob1]= distance_distribution(u1,mu1); 
    [u2_dist,prob2]= distance_distribution(u2,mu2);
    
   
    
    d = dW_Rp(u1_dist,u2_dist,prob1,prob2,p);
    
end