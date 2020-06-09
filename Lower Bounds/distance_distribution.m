function [d_dist,prob] = distance_distribution(d,mu)
% The function returns the distribution of distances of the metric measure space
% induced by the distance matrix d and the probability vector mu.

n= length(mu);
distance_probs = zeros(n);

for i = 1:n
    for j= 1:n
        distance_probs(i,j)=mu(i)*mu(j);
    end
end

v = sort(reshape(d,[],1));


d_dist = unique(v);
prob=zeros(length(d_dist),1);

for i = 1:length(d_dist)
	indices=find(d==d_dist(i));
	prob(i)=sum(distance_probs(indices));
end

    
    
end
