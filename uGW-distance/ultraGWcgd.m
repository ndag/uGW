function [result_final,coupl_final,result_vec,result_coupling_cell] = ultraGWcgd(ux,uy,mu_x,mu_y,p,iterations,num_samples, num_skips)
% This function approximates the ultrametric Gromov-Wasserstein distance
% via conditional gradient descent (Frank-Wolfe algorithm with fixed decreasing step size)

%   input: ux,uy      ultrametric distance matrices
%          mu_x,mu_y  weight (row)vectors 
%          p          order of the GW-distance 
%          iterations number of iterations 

%
%   output: 
%   
%           result          The ultrametric Gromov-Wasserstein distance approximated by
%                           conditional gradient descent
%           pi              The optimal coupling approximated by gradient
%                           descent

if num_samples >= 1
    [A,b]=gw_equality_constraints(mu_x,mu_y);
    coupls = coupling_ensemble(A,b,mu_x,mu_y,num_samples,num_skips);
    coupls=[ { mu_y'*mu_x};coupls];
else
    coupls = { mu_y'*mu_x};
end

result_vec = zeros(1,length(coupls));
result_coupling_cell=cell(length(coupls));

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
    result_coupling_cell{index}=pi;
    
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
[result_final, min_index]=min(result_vec);
coupl_final= result_coupling_cell{min_index};

end