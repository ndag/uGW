
load Us.mat

L = length(Us);
scale = [0 0.2 0.4 0.6];%Noise levels
S = length(scale);
    Yi = Us{35};
    uyi = Yi.ux;
    nyi = length(uyi);  
    muyi= ones(1,length(uyi))/length(uyi);

%parameters for uGWapprox
p= 1;
iterations= 10;
num_samples=5;
num_skips=20;
parpool(4)

for k = 1:S
    k
    uSLB= zeros(L,L);
    dSLB = zeros(L,L);
    dGW = zeros(L,L);
    uGW = zeros(L,L);
    lab = zeros(1,L); % labels


    for i=1:L
        Xi = Us{i};
        ui = Xi.ux;
        lab(i) = Xi.nb; % use number of blocks as label
        mui= ones(1,length(ui))/length(ui);
        ui = disturb_ultrametric_space(ui, scale(k));
        i
        for j=i+1:L
            Xj = Us{j};
            uj = Xj.ux;
            muj= ones(1,length(uj))/length(uj);
            uj = disturb_ultrametric_space(uj, scale(k));
            % calculate the lower bound/GW-dist etc.
            uSLB(i,j)= ugwslb(ui,uj,mui,muj,p);
            uGW(i,j)= ultraGWcgd(ui,uj,mui,muj,p,iterations,num_samples, num_skips);
            dSLB(i,j)= dgwslb(ui,uj,mui,muj,p);
            dGW(i,j)= dGWcgd(ui,uj,mui,muj,p,iterations,num_samples, num_skips);
        end
    end

    uSLB= max(uSLB,uSLB');
    uGW = max(uGW,uGW');
    dSLB= max(dSLB,dSLB');
    dGW = max(dGW,dGW');
    delete(gcp('nocreate'))

    noise = scale(k);
    
    [sl il] = sort(lab);
    subplot(4,S,k,'align'), imagesc(uSLB(il,il)); title(['uSLB perturbation at level ',num2str(noise)])
    axis square
    subplot(4,S,S+k,'align'), imagesc(uGW(il,il)); title(['uGW perturbation at level ',num2str(noise)])
    axis square
    subplot(4,S,2*S+k,'align'), imagesc(dSLB(il,il)); title(['dSLB perturbation at level ',num2str(noise)])
    axis square
    subplot(4,S,3*S+k,'align'), imagesc(dGW(il,il)); title(['dGW perturbation at level ',num2str(noise)])
    axis square

end

