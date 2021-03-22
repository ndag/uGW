
load Us_ind.mat

L = length(Us);

dSLB= zeros(L,L);
dTLB = zeros(L,L);
dGW = zeros(L,L);
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


        % Calculate the lower bounds/GW-dist etc.
        dSLB(i,j)= dgwslb(ui,uj,mui,muj,p);
        dGW(i,j)= dGWcgd(ui,uj,mui,muj,p,iterations,num_samples, num_skips);
        dTLB(i,j)= dgwtlb(ui,uj,mui,muj,p);
        
    end
end

dSLB= max(dSLB,dSLB');
dTLB = max(dTLB,dTLB');
dGW = max(dGW,dGW');

[~, il] = sort(lab);

dSLB=dSLB(il,il);
dTLB=dTLB(il,il);
dGW=dGW(il,il);
delete(gcp('nocreate'))