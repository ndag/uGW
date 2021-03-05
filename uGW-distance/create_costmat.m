
function cost_mat = create_costmat(ux,uy,coupling,p)
% The function calculates the cost matrix required for the gradient step in
% UltraGWapprox

n=length(ux);
m=length(uy);


cost_mat=zeros(n,m);




for i=1:n
    for j=1:m
        tmp=0;
        for k=1:n
            for l=1:m
                tmp=tmp+(delta_infinity(ux(i,k),uy(j,l)))^p*coupling(k,l);
            end
        end
        cost_mat(i,j)=2*tmp;
    end
end

end










































