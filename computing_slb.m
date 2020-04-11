load Us.mat

L = length(Us);

dm = zeros(L,L); %ugwslb
dm1 = zeros(L,L); %dgwslb
lab = zeros(1,L); % labels

for i=1:L
    Xi = Us{i};
    ui = Xi.ux;
    lab(i) = Xi.nb; % use number of blocks as label
    %diam1X = sum(ui(:))/length(ui)^2;
    i
    for j=i+1:L
        Xj = Us{j};
        uj = Xj.ux;
        %diam1Y = sum(uj(:))/length(uj)^2;
        
        % here calculate the lower bound (right now computing difference betwene 1-diams)
        
        %dm(i,j) = abs(diam1X-diam1Y);
        dm(i,j) = ugwslb(ui,uj,1);
        dm1(i,j) = dgwslb(ui,uj);
        
    end
end

dm = max(dm,dm');
dm1 = max(dm1,dm1');
[sl il] = sort(lab);
subplot(1,2,1,'align'), imagesc(dm(il,il)); title('uGW-slb')
subplot(1,2,2,'align'), imagesc(dm1(il,il)); title('dGW-slb')
axis square

