
load Us_ind.mat

L = length(Us);

uSLB= zeros(L,L);
uTLB = zeros(L,L);
uGW = zeros(L,L);
lab = zeros(1,L); % labels


%parameters for uGWapprox
p= 2;
iterations= 200;
num_samples=5;
num_skips=20;
parpool(4)
for i=1:L
    Xi = Us{i};
    ui = Xi.ux;
    lab(i) = Xi.nb; % use number of blocks as label
    mui= ones(1,length(ui))/length(ui);
    for j=i+1:L
        j
        Xj = Us{j};
        uj = Xj.ux;
        muj= ones(1,length(uj))/length(uj);


        % calculate the lower bound/GW-dist etc.
        uSLB(i,j)= ugwslb(ui,uj,mui,muj,p);
        uGW(i,j)= ultraGWcgd(ui,uj,mui,muj,p,iterations,num_samples, num_skips);
        uTLB(i,j)= ugwtlb(ui,uj,mui,muj,p);
        
    end
end

uSLB= max(uSLB,uSLB');
uTLB = max(uTLB,uTLB');
uGW = max(uGW,uGW');

[~, il] = sort(lab);

uSLB=uSLB(il,il);
uTLB=uTLB(il,il);
uGW=uGW(il,il);
delete(gcp('nocreate'))

