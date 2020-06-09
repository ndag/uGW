% The Function ultraWasserstein calculates the p-th power of the p-Wasserstein distance between two measures mu and nu (both vectors) on an ultrametric space (encoded by a distance matrix). The 
% implementation is relatively efficient and seems to perform well.  

% A first example

dist_mat = [0,1,4,4,4,4,4;1,0,4,4,4,4,4;4,4,0,0.5,4,4,4;4,4,0.5,0,4,4,4;4,4,4,4,0,3,3;4,4,4,4,3,0,3;4,4,4,4,3,3,0];
mu =  [0.1,0.1,0.2,0.05,0.05, 0.3,0.2];
nu = ones(1,7)/7;

mex ultraWasserstein.cpp -larmadillo

p=1;

result = ultraWasserstein(mu,nu,dist_mat,p)

% A big, randomly generated metric space

dist_mat_two  = importdata('ux5000.txt');

dist_mat_two=round(10*dist_mat_two,2);
% Attention! ux10000+diag(10000) contains zeros!
% This leads to errors! Thus:

dist_mat_two=dist_mat_two+(ones(5000,5000)-eye(5000));

mex ultraWasserstein.cpp -larmadillo

mu_two =  abs(normrnd(0,1,[1,5000]));
mu_two = mu_two/sum(mu_two);

nu_two =  abs(normrnd(0,1,[1,5000]));
nu_two = nu_two/sum(nu_two);
p=1;

tic
result = ultraWasserstein(mu_two,nu_two,dist_mat_two,p)
toc
% Comparison to the usage of the linkage function
tic
linkage(dist_mat_two);
toc

% On our server the ultraWasserstein function is 80 times faster than the 
% standard linkage function.
load('dist_mat_three.mat')
mu_three =  abs(normrnd(3,10,[1,5000]));
mu_three = mu_three/sum(mu_three);

nu_three =  abs(normrnd(3,10,[1,5000]));
nu_three = nu_three/sum(nu_three);
p=1;

tic
result = ultraWasserstein(mu_three,nu_three,dist_mat_three,p)
toc
% Comparison to the usage of the linkage function
tic
linkage(dist_mat_three);
toc

% The 1-D case
num = 5

pos = sort(abs(normrnd(3,10,[1,num])));

mu_four =  abs(normrnd(3,10,[1,num]));
mu_four = mu_four/sum(mu_four);

nu_four =  abs(normrnd(3,10,[1,num]));
nu_four = nu_four/sum(nu_four);

dist_mat_four = zeros(num);
for i = 1:num
    for j = 1:num
        if(pos(i) ~=pos(j))
            dist_mat_four(i,j)=max(pos(i),pos(j));
        end
    end
end
result = (ultraWasserstein(mu_four,nu_four,dist_mat_four,1))
[dist_out,gamma_out] = mexEMD(nu_four.',mu_four.',dist_mat_four)
uW_R(sort(pos),mu_four,nu_four,1)
