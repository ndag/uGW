# uGW

## ugwslb.m 
**ugwslb.m** is a Matlab function computing the second lower bound of uGW of two ultrametric measure spaces.

**Syntax:**
```
% d = ugwslb(u1,u2,mu1,mu2,p)
%   u1  - ultrametric distance matrix
%   u2  - ultrametric distance matrix
%   mu1 - probability vector
%   mu2 - probability vector
%   p   - real number >=1
%
% Returns:
%   d - pth-uGW-second lower bound between u1 and u2 with probability measures mu1 and mu2, respectively.
%   plan  - the corresponding transport plan
```
## ugwtlb.m 
**ugwtlb.m** is a Matlab function computing the third lower bound of uGW of two ultrametric measure spaces.

**Syntax:**
```
% [res,coupl] = ugwtlb(u1,u2,mu1,mu2,p)
%   u1  - ultrametric distance matrix
%   u2  - ultrametric distance matrix
%   mu1 - probability vector
%   mu2 - probability vector
%   p   - real number >=1
%
% Returns:
%   res    - pth-uGW-third lower bound between u1 and u2 with probability measures mu1 and mu2, respectively.
%   coupl  - the corresponding optimal coupling
```
## ultraGWcgd.m 
**ultraGWcgd.m** is a Matlab function that approximates the ultrametric Gromov-Wasserstein distance via conditional gradient descent (Frank-Wolfe algorithm with fixed decreasing step size). It is possible to choose multtiple random starting points for the gradient descent algorithm.

**Syntax:**
```
% [res,coupl,result_vec,result_coupl_cell] = ultraGWcgd(ux,uy,mux,muy,p,iterations,num_samples, num_skips)
%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%   iterations - number of iterations for the gradient descent
%   num_samples - number of random starting points (if this is set to 0 only the independece coupling is used)
%   num_skips - number of skips between two random couplings 
%   
% Returns:
%   res   - the best approximation of the ultrametric Gromov-Wasserstein distance of order p between u1 and u2 with probability measures mu1 and mu2, respectively.
%   plan  - the corresponding coupling
%   result_vec - a vector containing all stationary points of the gradient descent started from the different couplings
%   result_plan_cell - the collection of corresponding couplings
```
## uW_R.m.m 
**uW_R.m** is a Matlab function that implements the Wasserstein distance on (R,Delta_infinity).

**Syntax:**
```
%  d = uW_R(pos,mu1,mu2,p)
%   mu1 - probability vector
%   mu2 - probability vector
%   p   - real number >=1

%   
% Returns:
%   d   - the (R,Delta_infinity)-Wasserstein distance between the probability measures mu1 and mu2.

```

## computing_slb.m 
**computing_slb.m** is a Matlab code computing uGW-slb and dGW-slb with p=1 between each pair of ultrametric spaces with uniform probability measures from a given file **Us.mat** and generates corresponding distance matrices sorted by number of blocks in ultrametric spaces. Before running the code, one should first run either one of the following Matlab codes to generate random ultrametric spaces recorded in file **Us.mat**:

* **test_subsample.m** Subsample ultrametric spaces from some given randomly generated ultrametric spaces.
* **test_independent.m** Generates ultrametric spaces independently.

## visual.m 
**visual.m** computes uGW-slb and dGW-slb with p=1 between any given pair (i,j) of spaces in **Us.mat** and generates the corresponding dendrograms and distance distributions.

**Syntax:**
```
% [d,d1] = visual(i,j)
%   i - index in Us.mat
%   j - index in Us.mat
%
% Returns:
%   d - uGW-second lower bound between Us{i} and Us{j}
%   d1 - dGW-second lower bound between Us{i} and Us{j}
```
