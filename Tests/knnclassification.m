addpath(genpath('uGW-Project'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% knn-classification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% k-nn classification
labels=[repmat(0,1,200),repmat(1,1,200)];

reps = 5000;
N=50;
k = 5;
num_miss_spezi=zeros(1,reps);

choice = 2;

if (choice ==1)
    load('uSLBInfluenzaTrees.mat')
    dist_mat = uSLB;

elseif(choice==2)
    load('SLBInfluenzaTrees.mat')
    dist_mat = SLB;
elseif(choice ==3)
    load('colijn_distmat.mat')
    dist_mat = colijn_distmat;
end


for i = 1:reps
    sample1 = datasample(1:200,N,'Replace',false);
    sample2 = datasample(201:400,N,'Replace',false);
    
    tmp1=dist_mat(sample1,:);
    tmp2=dist_mat(sample2,:);
    
    tmp1(:,[sample1,sample2])=[];
    tmp2(:,[sample1, sample2])=[];
    
    new_dist_mat = [tmp1;tmp2];
    drawn_labels= labels([sample1,sample2]);
    
    dim=size(new_dist_mat);
    rem_cols = 1:400;
    rem_cols([sample1,sample2])=[];
    
    for j = 1:(dim(2))
        [tmp3,tmp_index] = sort(new_dist_mat(:,j));

        new_label = 0;
        if(nnz(drawn_labels(tmp_index(1:k)))>k/2)
            new_label=1;
        end
        if(new_label ~= labels(rem_cols(j)))
            num_miss_spezi(i)= num_miss_spezi(i)+1;
        end
    end
 
end


mean(num_miss_spezi)




































