
load Us.mat

L = length(Us);
scale = [0 0.2 0.4 0.6];%Noise levels
S = length(scale);


DistUs = cell(L,S);

% Preparation of the disturbed spaces
for k = 1:S
     for i=1:L
        Y.uy = disturb_ultrametric_space(Us{i}.ux, scale(k));
        Y.muy =  ones(1,length(Us{i}.ux))/length(Us{i}.ux);
        DistUs{i,k} = Y;
     end

end

save DisturbedUs.mat DistUs


%parameters for uGWapprox
p= 1;
iterations= 100;
num_samples=5;
num_skips=20;
%parpool(4)

 uSLB_res = cell(S);
 dSLB_res = cell(S);
 uTLB_res = cell(S);
 dTLB_res = cell(S);
 dGW_res  = cell(S);
 uGW_res  = cell(S);


for k = 1:S
    k
    uSLB= zeros(L,L);
    dSLB = zeros(L,L);
    uTLB= zeros(L,L);
    dTLB = zeros(L,L);
    dGW = zeros(L,L);
    uGW = zeros(L,L);


    for i=1:L
        i
        ui  = DistUs{i,k}.uy;
        mui = DistUs{i,k}.muy;

        for j=i:L

            uj  = DistUs{j,k}.uy;
            muj = DistUs{j,k}.muy;
            % calculate the lower bound/GW-dist etc.
            uSLB(i,j)= ugwslb(ui,uj,mui,muj,p);
            uTLB(i,j)= ugwtlb(ui,uj,mui,muj,p);
            uGW(i,j)= ultraGWcgd(ui,uj,mui,muj,p,iterations,num_samples, num_skips);
            dSLB(i,j)= dgwslb(ui,uj,mui,muj,p);
            dTLB(i,j)= dgwtlb(ui,uj,mui,muj,p);
            dGW(i,j)= dGWcgd(ui,uj,mui,muj,p,iterations,num_samples, num_skips);
        end
    end

    uSLB= max(uSLB,uSLB');
    uGW = max(uGW,uGW');
    dSLB= max(dSLB,dSLB');
    dGW = max(dGW,dGW');
    dTLB = max(dTLB,dTLB');
    uTLB = max(uTLB,uTLB');

    %delete(gcp('nocreate'))
    
 uSLB_res{k} = uSLB;
 dSLB_res{k} = dSLB;
 uTLB_res{k} = uTLB;
 dTLB_res{k} = dTLB;
 dGW_res{k}  = dGW;
 uGW_res{k}  = uGW;

end

save('uGWApprox','uSLB_res','dSLB_res', 'uTLB_res','dTLB_res', 'uGW_res','dGW_res')

