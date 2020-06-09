function [d_dist,prob] = distance_distribution_fast(d,mu)
% The function returns the distribution of distances of the metric measure space
% induced by the distance matrix d and the probability vector mu.

n= length(mu);
distance_probs = zeros(n);

for i = 1:n
    for j= 1:n
        distance_probs(i,j)=mu(i)*mu(j);
    end
end


[dist_vec,dist_vec_order] = sort(reshape(d,[],1));
tmp = reshape(distance_probs,[],1);
prob_vec = tmp(dist_vec_order);

num_of_unique_el =n*(n-1)/2+1;
d_dist_tmp = zeros(num_of_unique_el ,1);
prob_tmp=zeros(num_of_unique_el ,1);

d_dist_tmp(1)=0;
prob_tmp(1)=sum(prob_vec(1:n));

j=n+1;
i= 2;
while j <= n^2
    d_dist_tmp(i) = dist_vec(j);
    prob_tmp(i)= prob_vec(j);
    k = j+1;
    while  k<=n^2 && d_dist_tmp(i) ==  dist_vec(k)
        prob_tmp(i)= prob_tmp(i)+prob_vec(k);
        k = k+1;
    end     
    j=k;
    i=i+1;
        
end

    d_dist=d_dist_tmp(1:(i-1));
    prob= prob_tmp(1:(i-1));
    
end
