function d = uGWinfty(ux, uy, mux, muy)

%   This function computes the ultrametric Gromov-Wasserstein distance of
%   order p=infinty.

%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%
% Returns:
%   res - infinity-uGW distance between ux and uy with probability measures mux and muy, respectively.

d=0;

spec = unique([unique(ux(:)); unique(uy(:))]);
spec = sort(spec,'descend');

ns = length(spec);

for c=1:ns
    t = spec(c);
    
    [subux,submux] = quotientUMFPS(ux,mux,t);
    [subuy,submuy] = quotientUMFPS(uy,muy,t);
    
    if ~is_iso(subux,subuy,submux,submuy)
        d = spec(max([c-1,1]));
        %d = spec(c);
        return
    end
end
