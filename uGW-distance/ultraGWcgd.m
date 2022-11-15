function [res,plan,result_vec,result_plan_cell] = ultraGWcgd(ux, uy, mu_x, mu_y, p, iterations, num_samples, num_skips)

% This function approximates the ultrametric Gromov-Wasserstein distance of order p via conditional gradient descent 
% (Frank-Wolfe algorithm with fixed decreasing step size). It is possible to start the gradient descent algorithm from multiple random couplings.

%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%   iterations - number of iterations for the gradient descent
%   num_samples - number of random starting couplings for the gradient descent (if this is set to 0 only the independece coupling is used)
%   num_skips - number of skips between two random couplings 
%   
% Returns:
%   res   - the best approximation of the ultrametric Gromov-Wasserstein distance of order p between ux and uy with probability measures mux and muy, respectively.
%   plan  - the corresponding coupling
%   result_vec - a vector containing all stationary points of the gradient descent started from the different couplings
%   result_plan_cell - the collection of corresponding couplings

if num_samples >= 1
    [A,b]=gw_equality_constraints(mu_x,mu_y);
    coupls = coupling_ensemble(A,b,mu_x,mu_y,num_samples,num_skips);
    coupls=[ { mu_y'*mu_x};coupls];
else
    coupls = { mu_y'*mu_x};
end

result_vec = zeros(1,length(coupls));
result_plan_cell=cell(1,length(coupls));

for index=1:length(coupls)
     pi_old = coupls{index};
     pi_old = pi_old';
    for i=1:iterations
    
        %cost=create_costmat(ux,uy,pi_old,p);
        if(p==1)
            cost=construct_cost_one(ux,uy,pi_old);
        else
            cost=construct_ultra_cost_mat(ux,uy,pi_old,p);
        end
        [~,ot_plan]=mexEMD(mu_x',mu_y',cost);
        step_size=2/(i+2);
        pi_new = (1-step_size)*pi_old+step_size*ot_plan;
        pi_old =pi_new;
    end
    pi= pi_new;
    result_plan_cell{index}=pi;
    
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
result_vec(index)=result^(1/p);
end
[res, min_index]=min(result_vec);
plan= result_plan_cell{min_index};

end