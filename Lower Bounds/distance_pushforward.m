function [d_dist,prob] = distance_pushforward(d,mu)
% The function returns the distribution of distances of the metric measure space
% induced by the distance matrix d and the probability vector mu.

[d_dist,~,tmp3]=unique(d);
prob = zeros(1,length(d_dist));
for i=1:length(mu)
    prob(tmp3(i))=prob(tmp3(i))+mu(i);
end

end