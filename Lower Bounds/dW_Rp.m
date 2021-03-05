function cost=dW_Rp(u, v,u_weights, v_weights, p)
    %Solves the Earth Movers distance problem between 1d measures and
    %returns the associated cost

    [u,tmp1] =sort(u);
    u_weights=u_weights(tmp1);

    [v,tmp2] =sort(v);
    v_weights=v_weights(tmp2);
    
    cost = 0;
    
    n = length(u_weights);
    m = length(v_weights);

    i = 1;
    w_i = u_weights(1);
    j = 1;
    w_j = v_weights(1);

    m_ij = 0;
    
    while(true)
        m_ij = abs(u(i) - v(j))^p;
        if (w_i < w_j || j == m)
            cost = cost+m_ij * w_i;
            i = i+1;
            if (i == n+1)
                break;
            end    
            w_j = w_j-w_i;
            w_i = u_weights(i);
        else
            cost = cost+ m_ij * w_j;
            j = j+1;
            if (j == m+1)
                break;
            end
            w_i =w_i- w_j;
            w_j = v_weights(j);
        end    
    end
   end