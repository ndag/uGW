function [projected_mu]=project_mu(mu,A,b,P,product_mu)

    % Input: coupling-shaped matrix mu; equality constraints A,b from gw_equality_constraints
    %        function; product coupling of some probability measures p and q.
    %        P is a projection matrix onto row space of A.
    % Output: Orthogonal projection of mu onto the affine subspace determined by A,b.

    [m,n]=size(product_mu);

    % Create the vector to actually project and reshape
    vec_to_project = mu - product_mu;
    vec_to_project = reshape(vec_to_project,n*m,1);

    % Project it
    vec_to_project = vec_to_project - P*vec_to_project;

    projected_mu = reshape(product_mu,m*n,1) + vec_to_project;

    projected_mu = reshape(projected_mu,m,n);

    end