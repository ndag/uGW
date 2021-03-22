function [res,plan] = dgwtlb(u1,u2,mu1,mu2, p)
% the function calculates the third lower bound of dgw between two
% ultrametric spaces with arbitrary measures
    m = length(mu1);
    n= length(mu2);
    costmat = zeros(m,n);
    
    u1_dist=cell(m);
    prob1=cell(m);
    
    u2_dist=cell(n);
    prob2=cell(n);
    
    for i =1:m
        [u1_dist{i},prob1{i}]= distance_pushforward(u1(i,:),mu1); 
    end
    
    for j = 1:n
        [u2_dist{j},prob2{j}]= distance_pushforward(u2(j,:),mu2);
    end
    
    for i =1:m
        for j = 1:n
            costmat(i,j)= dW_Rp(u1_dist{i},u2_dist{j},prob1{i},prob2{j},p);
        end
    end

    
    [res,plan] = mexEMD(mu1.',mu2.',costmat);
    res = res^(1/p);
end
