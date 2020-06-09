function d = ugwslb_general_old(u1,u2,mu1,mu2, p)
%the function calculates the second lower bound of ugw between two
%ultrametric spaces with arbitrary measures
    
    [u1_dist,prob1]= distance_distribution(u1,mu1); 
    [u2_dist,prob2]= distance_distribution(u2,mu2);
    
    %generate position by taking union
    pos = union(u1_dist,u2_dist);
    pos = pos';
    

    n = length(pos);
    w1 = zeros(1,n);
    w2 = zeros(1,n);

    for i = 1: length(prob1)
        w1(pos == u1_dist(i)) = prob1(i);
    end

    for i = 1: length(prob2)
        w2(pos == u2_dist(i)) = prob2(i);
    end
    
    
    d = uW_R(pos,w1,w2,p);
    
end