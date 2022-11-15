function [res,plan] = ugwtlb(ux, uy, mux, muy, p)

%   A function for computing the third lower bound of uGW of two ultrametric measure spaces
%
%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res    - pth-uGW-third lower bound between ux and uy with probability measures mux and muy, respectively.
%   coupl  - the corresponding optimal coupling

    m = length(mux);
    n= length(muy);
    costmat = zeros(m,n);

    u1_dist=cell(m);
    prob1=cell(m);
    
    u2_dist=cell(n);
    prob2=cell(n);
    
    for i =1:m
        [u1_dist{i},prob1{i}]= distance_pushforward(ux(i,:),mux); 
    end
    
    for j = 1:n
        [u2_dist{j},prob2{j}]= distance_pushforward(uy(j,:),muy);
    end
    
    for i =1:m
        for j = 1:n
            costmat(i,j)=utlb_cost_mat_entry(u1_dist{i},prob1{i},u2_dist{j},prob2{j},p);
        end
    end
    
    
    [res,plan] = mexEMD(mux.',muy.',costmat);
    res =(res)^(1/p);
end



