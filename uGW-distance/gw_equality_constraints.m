function [A,b] = gw_equality_constraints(mux,muy)
% Inputs: probability row vectors
% Output: matrices A and b defining equality constraints

    m = length(mux);
    n = length(muy);

    A_p_type = kron(eye(m),ones(1,n));


    A_q_type=kron(ones(1,m),eye(n));
    
    A = vertcat(A_p_type,A_q_type);
    b =  horzcat(mux,muy);
    
end