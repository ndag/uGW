function res = uGWinfty(ux,uy,mux,muy)

res=0;

spec = unique([unique(ux(:)); unique(uy(:))]);
spec = sort(spec,'descend');

ns = length(spec);

for c=1:ns
    t = spec(c);
    
    [subux,submux] = quotientUMFPS(ux,mux,t);
    [subuy,submuy] = quotientUMFPS(uy,muy,t);
    
    if ~is_iso(subux,subuy,submux,submuy)
        res = spec(max([c-1,1]));
        %d = spec(c);
        return
    end
end
