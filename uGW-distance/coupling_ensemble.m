function [Markov_steps]=coupling_ensemble(A,b,mux,muy,num_samples,num_skips,mu_initial)

    % Inputs: equality constraints A,b; probability vectors p,q; number of steps
    %         to take in the Markov chain; initialization
    % Output: Ensemble of couplings from the probability simplex.

   product_mu = muy'*mux;

    if (nargin<7)
        mu_initial = product_mu;
    end

    % Find orthonormal basis for row space of A
    Q = orth(A');
    % Create projector onto the row space of A
    P = Q*Q';

    num_steps = num_samples*num_skips;

    Markov_steps = cell(num_samples,1);

    for j = 1:num_steps
        mu_new = markov_hit_and_run_step(A,b,P,mux,muy,mu_initial);
        mu_initial = mu_new;
        if (mod(j,num_skips)==0)
            Markov_steps{j/num_skips}=mu_new;
        end
    end
    
end