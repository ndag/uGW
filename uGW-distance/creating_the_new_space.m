function [mu1,mu2, dist_mat] = creating_the_new_space(ux,uy, mu_x,mu_y, points_to_be_id)
% This function creates the common metric spaces required to implement the
% GW-representation from Thm. 4.3
%
%   input: ultrametric distance matrices ux,uy
%          weight vectors mu_x, mu_y
%          l x 2 matrix of the indices of the l points to be identified
%          Ex.: Consider the spces {x_1,x_2,x_3} and {y_1,y_2}
%          In order to identify x_1 and y_2 points_to_be_id=[1,2]
%
%   output: two measures on the common space mu1,mu2
%           distance matrix of the common space    
m = length(mu_x);
n=length(mu_y);
num_of_points_tbid = size(points_to_be_id,1);
k = n+m-num_of_points_tbid;

mu1= zeros(k,1);
mu2=zeros(k,1);
dist_mat = zeros(k,k);

ptr_y = points_to_be_id(:,2);
rem_ind_y= setdiff(1:n,ptr_y);

mu1(1:m)=mu_x;
mu2(points_to_be_id(:,1))=mu_y(points_to_be_id(:,2));
mu2(m+1:k)=mu_y(rem_ind_y);

dist_mat(1:m,1:m)=ux;
dist_mat((m+1):k,(m+1):k)=uy(rem_ind_y,rem_ind_y);


for i = 1:num_of_points_tbid
    dist_mat(points_to_be_id(i,1),(m+1):k)=uy(points_to_be_id(i,2),rem_ind_y);
    dist_mat((m+1):k,points_to_be_id(i,1))=uy(points_to_be_id(i,2),rem_ind_y);
end


for j=setdiff(1:m,points_to_be_id(:,1))
    counter = 1;
    for h = setdiff(1:n,points_to_be_id(:,2))
        dist_mat(j,m+counter)= min(max(ux(j,points_to_be_id(:,1)),uy(h,points_to_be_id(:,2)))); 
        dist_mat(m+counter,j)=dist_mat(j,m+counter);
        counter = counter +1;
    end
end
end

