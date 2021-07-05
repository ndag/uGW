function disturbed_um_space = disturb_ultrametric_space(ux, scale)

disturbed_um_space = ux;
if scale ~= 0
    [cluster, min_dist]= find_cluster(ux, scale);
        for j = 1:length(min_dist)
            if(length(cluster{j})>1)
                ux_tmp = disturb(ux(cluster{j},cluster{j}), scale);
                disturbed_um_space(cluster{j},cluster{j})= ux_tmp;
            end
        end
end
end


function [cluster_elem, min_dist] = find_cluster(ux, t)

tmp_clus = squareform(ux);
link = linkage(tmp_clus);
clusters = cluster(link,'Criterion','distance','Cutoff',t);

noc = max(clusters);

cluster_elem = cell(1,noc);

represent = zeros(noc,1);

for i = 1:noc
    tmp = find(clusters == i);
    cluster_elem{i}=tmp;
    represent(i)=tmp(1);
end

if i>1
    tmp_dist = ux(represent,represent);
    tmp_dist=tmp_dist+diag(Inf(size(diag(tmp_dist))));
    min_dist = min(tmp_dist);
else
    min_dist = 0;
end

end


function ux_disturbed = disturb(ux, scale)
r = 1;%ratio of noise
ux_disturbed=ux;
[spectrum,~,r2] = unique(ux);

error = sort(unifrnd(0,abs(scale-spectrum(end)),length(spectrum),1));

for j = 1:length(r2)
    if(r2(j)>1)
        ux_disturbed(j)=ux_disturbed(j)+r * error(r2(j));
    end
end

end