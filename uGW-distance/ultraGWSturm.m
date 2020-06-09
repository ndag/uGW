function [result,com_space_dist, id_points] = ultraGWSturm(ux,uy,p)
% This function calculates the ultrametric Gromov-Wasserstein distance
% between two ultrametric spaces equipped with the uniform measure
% by iterating over all possible spaces.  This is extremely inefficient!!! 
% It is computationally only feasible for small spaces!
%
%   input: ux,uy    ultrametric distance matrices 
%          p        order 
%
%   output: 
%   
%   result          The ultrametric GW-distance between (ux,mu_x) and (uy,mu_y) 
%   com_space_dist  The distance matrix of the optimal common space
%   id_points       l x 2 matrix of the indices of the l points to be identified
%                   in the optimal common space
%                   Ex.: Consider the spaces {x_1,x_2,x_3} and {y_1,y_2}
%                   If x_1 and y_2 are identified then id_points=[1,2]

% Remarks: 1. The calculation becomes even more messy for general measures.
%          (see comment in the code).

m = length(ux);
n= length(uy);
result = Inf;
mu_x=ones(m,1)/m;
mu_y=ones(n,1)/n;

for i=1:m
    for j=1:n
        [tmp,tmp_dist]= ultraGWSturmaux(ux,uy,mu_x,mu_y,[i,j],p);
        if(tmp<result)
            result = tmp;
            com_space_dist = tmp_dist;
            id_points=[i,j];
        end
    end
end

for k=2:min(n,m)
    comb_x= nchoosek(1:m,k);
    comb_y= nchoosek(1:n,k);
    s_x =size(comb_x);
    s_y=size(comb_y);
    for i= 1:(s_x(1))
        for j = 1:(s_y(1))
            % This is the point, where things get even more troublesome.
            % It is not sufficent to identify one isomorphism (which could be done by Find_Corr or an
            % adapted version of uGH), but one as to further iterate over all possible isomorphisms.
            
            d=Find_Corr(ux(comb_x(i,:),comb_x(i,:)),uy(comb_y(j,:),comb_y(j,:)),10^(-10)); % Here, it might be possible to
                                                                                           % speed up the code by adapting 
                                                                                           % the function uGH. 
            if(length(d(:,1))>1)
                [rId, cId] = find(d);
                rId= comb_x(i,rId);
                cId= comb_y(j,cId);
                [tmp,tmp_dist]= ultraGWSturmaux(ux,uy,mu_x,mu_y,[rId',cId'],p);
                 if(tmp<result)
                    result = tmp;
                    com_space_dist = tmp_dist;
                    id_points=[rId',cId'];
                end
            end
        end
    end
end
end

