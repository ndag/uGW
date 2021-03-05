function [result_final,coupl_final,result_vec,result_coupling_cell] = dGWcgd(dx,dy,mu_x,mu_y,p,iterations,num_samples, num_skips)
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
    
        %cost=create_costmat_GW(dx,dy,pi_old,p);
        if(p==1)
            cost=construct_cost_one_dGW(dx,dy,pi_old);
        else
            cost=construct_cost_mat_dGW(dx,dy,pi_old,p);
        end
        [~,ot_plan]=mexEMD(mu_x',mu_y',cost);
        step_size=2/(i+2);
        pi_new = (1-step_size)*pi_old+step_size*ot_plan;
        pi_old =pi_new;
    end
    pi= pi_new;
    result_coupling_cell{index}=pi;
    
    result= 0;
    for i =1:length(dx)
        for j=1:length(dy)
            for k=1:length(dx)
                for l=1:length(dy)
                    result= result + pi(i,j)*pi(k,l)*(abs(dx(i,k)-dy(j,l)))^p;
                end     
            end
        end
    end
result_vec(index)=result^(1/p);
end
[result_final, min_index]=min(result_vec);
coupl_final= result_coupling_cell{min_index};

end
