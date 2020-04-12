function [d,d1] = visual(i,j)
    %dm is the slb matrix
    load Us.mat;
    Xi = Us{i};
    ui = Xi.ux;
    nxi = length(ui);
    Xj = Us{j};
    uj = Xj.ux;
    nxj = length(uj);
    
    d = ugwslb(ui,uj,1);
    d1 = dgwslb(ui,uj);
    
    lnkxi = linkage(squareform(ui));
    subplot(2,2,1), dendrogram(lnkxi); title(['Xi ugw-slb=' num2str(d)]);
    axis([0 nxi+1 0 0.6])
    %axis([0 ny 0 diam])



    lnkxj = linkage(squareform(uj));
    subplot(2,2,2), dendrogram(lnkxj); title(['Xj dgw-slb=' num2str(d1) ]);
    axis([0 nxj+1 0 0.6])
    %axis([0 nx 0 diam])
    
    subplot(2,2,3), 
    hi = histogram(squareform(ui),'Normalization','probability');
    hi.NumBins = 100;
    subplot(2,2,4), 
    hj = histogram(squareform(uj),'Normalization','probability');
    hj.NumBins = 100;
end