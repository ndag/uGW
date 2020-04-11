function d = visual(i,j)
    %dm is the slb matrix
    load Us.mat;
    Xi = Us{i};
    nbi = Xi.nb
    UXi = Xi.ux;
    ui = Xi.ux;
    Xj = Us{j};
    nbj = Xj.nb
    UXj = Xj.ux;
    uj = Xj.ux;
    
    d = ugwslb(ui,uj,1);
    d1 = dgwslb(ui,uj);
    
    lnkxi = linkage(squareform(UXi));
    subplot(2,2,1), dendrogram(lnkxi); title(['Xi ugw-slb=' num2str(d)])
    %axis([0 ny 0 diam])



    lnkxj = linkage(squareform(UXj));
    subplot(2,2,2), dendrogram(lnkxj); title(['Xj dgw-slb=' num2str(d1) ]);
    %axis([0 nx 0 diam])
    
    subplot(2,2,3), 
    hi = histogram(squareform(UXi),'Normalization','probability');
    hi.NumBins = 100;
    subplot(2,2,4), 
    hj = histogram(squareform(UXj),'Normalization','probability');
    hj.NumBins = 100;
end