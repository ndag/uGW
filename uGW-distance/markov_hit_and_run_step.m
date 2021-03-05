function [mu_new]=markov_hit_and_run_step(A,b,P,mux,muy,mu_initial)
    % Input: equality constraints A,b from gw_equality_constraints; pair of
    %       probability vectors p, q; initialization
    %        P is a projection matrix onto row space of A.
    % Output: new coupling measure after a hit-and-run step.

    m = length(mux);
    n = length(muy);

    product_mu = muy'*mux;

    if (nargin<6)
        mu_initial = product_mu;
    end
    
    mu_initial = project_mu(mu_initial,A,b,P,product_mu);
    % Project to the affine subspace
    % We assume mu_initial already lives there, but this will help with accumulation of numerical error

    mu_initial = reshape(mu_initial,m*n,1);

    % Choose a random direction
    direction = normrnd(0,1, m*n,1);

    % Project to subspace of admissible directions

    direction = direction - P*direction;

    % Renormalize

    direction = direction/norm(direction,'fro');

    % Determine how far to move while staying in the polytope - These are inequality bounds,
    % so we just need the entries to stay positive

    pos = direction > 1e-6;
    neg = direction < -1e-6;

    direction_pos = direction(pos);
    direction_neg = direction(neg);
    mu_initial_pos = mu_initial(pos);
    mu_initial_neg = mu_initial(neg);

    lower = max(-mu_initial_pos./direction_pos);
    upper = min(-mu_initial_neg./direction_neg);

    % Choose a random distance to move
    r = (upper - lower)*rand() + lower;

    mu_new = mu_initial + r*direction;
    mu_new = reshape(mu_new,n,m);

    end