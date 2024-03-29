

% NBs, nb: number of blocks
% NCs, nc: number of points per block
% SEPs, sep: minimal separatation between different blocks
% nb and nc are used to generate a randome ultrametric space Z with nb
% blocks and nc points per block
% X is generated from Z by subsampling points:
% NXs, nx: number of points from Z that will be sampled to form X

% results wil be saved into a cell array (see below for 'filename')

NXs = [10 10 10 20 20 20 30 30 30]
%NXs = [30 10 10 20 20 20 30 30 30]
NBs = [1 2 3 4 5]
NCs = [100]
SEPs = [1.5]
filename = 'Us.mat'

n1 = length(NXs);
n2 = length(NBs);
n3 = length(NCs);
n4 = length(SEPs);


N = n1*n2*n3*n4;

Us = cell(N,1);

for i=1:n1*n2*n3*n4
    [ix ib ic isep] = ind2sub([n1 n2 n3 n4], i);
    
    nb = NBs(ib);
    nc = NCs(ic);
    sep = SEPs(isep);
    nx = NXs(ix);
    
    Z = [];
    for k=1:nb
        Z = [Z; (k-1)*sep + rand(nc,1)];
    end
    
    nz = length(Z);
    dz = L2_distance(Z',Z');
    uz = max_subdominant_ultra(dz);
    rx = randperm(nz);
    Ix = rx(1:nx);
    ux = uz(Ix,Ix);
    ux = ux.*(1-eye(length(ux)));
    
    X.ux = ux;
    X.sep = sep;
    X.nb = nb;
    X.nc = nc;
    
    Us{i} = X;
    i
    
    %lnkx = linkage(squareform(ux));
    %dendrogram(lnkx,0); title('X')
    %axis([0 nx 0 max(ux(:))])
    %title(sprintf('nx = %d, nb = %d nc = %d sep = %d', nx, nb, nc, sep))
    %drawnow
end

save Us.mat Us



