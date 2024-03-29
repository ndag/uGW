function res = utlb_cost_mat_entry(u1_dist,prob1,u2_dist,prob2,p)
%the function calculates one cost matix entry for the calculation ofthe third lower bound of ugw between two
%ultrametric spaces with arbitrary measures
    
    %generate position by taking union
    pos = union(u1_dist,u2_dist);
    

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
    
    res = uW_R(pos,w1,w2,p);
    
end