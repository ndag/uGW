function [result,com_space_dist] = ultraGWSturmaux(ux,uy,mu_x,mu_y,points_to_be_id,p)
% This function calculates the Wasserstein distance the on th common metric
% space given by points_to_be_id (see the ultraGW-representation from Thm. 
% 4.3)
%   input: ultrametric distance matrices ux,uy
%          weight vectors mu_x, mu_y
%          l x 2 matrix of the indices of the l points to be identified
%          Ex.: Consider the spces {x_1,x_2,x_3} and {y_1,y_2}
%          In order to identify x_1 and y_2 points_to_be_id=[1,2]
%
%   output: 
%   
%   result          The Wasserstein distance between the measures mu_x and mu_y on 
%                   the common ultrametric space induced by points_to_be_id
%
%   com_space_dist  The distance matrix of the common space induced by
%                   points_to_be_id

[mu,nu,com_space_dist]=creating_the_new_space(ux,uy,mu_x,mu_y, points_to_be_id);
%mex ultraWasserstein.cpp -larmadillo; (if not executed somewhere else) 
result = (ultraWasserstein(mu,nu,com_space_dist,p))^(1/p);
end

