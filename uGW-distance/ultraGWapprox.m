function [result,pi] = ultraGWapprox(ux,uy,mu_x,mu_y,p,iterations)
% This function approximates the ultrametric Gromov-Wasserstein distance
% via conditional gradient descent (Frank-Wolfe algorithm with fixed decreasing step size)

%   input: ux,uy      ultrametric distance matrices
%          mu_x,mu_y  weight (row)vectors 
%          p          order of the GW-distance 
%          iterations number of iterations 

%
%   output: 
%   
%           result          The Gromov-Wasserstein distance approximated by
%                           conditional gradient descent
%           pi              The optimal coupling approximated by gradient
%                           descent

pi_old = mu_x'*mu_y;

for i=1:iterations
    
    cost=create_costmat(ux,uy,pi_old,p);
    [~,ot_plan]=mexEMD(mu_x',mu_y',cost);
    step_size=2/(i+2);
    pi_new = (1-step_size)*pi_old+step_size*ot_plan;
    pi_old =pi_new;
end


pi= pi_new;
result= 0;

for i =1:length(ux)
    for j=1:length(uy)
        for k=1:length(ux)
            for l=1:length(uy)
                result= result + pi(i,j)*pi(k,l)*(delta_infinity(ux(i,k),uy(j,l)))^p;
            end     
        end
    end
end
result=result^(1/p);
end