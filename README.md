This is an implementation of notions related to the ultrametric Gromov-Wasserstein distance developed by [Facundo Memoli](https://people.math.osu.edu/memolitechera.1/), [Axel Munk](http://www.stochastik.math.uni-goettingen.de/munk), [Zhengchao Wan](https://math.osu.edu/people/wan.252-0) and Christoph Weitkamp in the paper https://arxiv.org/abs/2101.05756. 

# uGW
## Lower Bounds
### ugwslb.m 
**ugwslb.m** is a Matlab function for computing the second lower bound of uGW of two ultrametric measure spaces.

**Syntax:**
```
% res = ugwslb(ux,uy,mux,muy,p)
%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res - pth-uGW-second lower bound between ux and uy with probability measures mux and muy, respectively.
```
### ugwtlb.m 
**ugwtlb.m** is a Matlab function for computing the third lower bound of uGW of two ultrametric measure spaces.

**Syntax:**
```
% [res,coupl] = ugwtlb(ux,uy,mux,muy,p)
%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res    - pth-uGW-third lower bound between ux and uy with probability measures mux and muy, respectively.
%   coupl  - the corresponding optimal coupling
```
### uW_R.m
**uW_R.m** is a Matlab function that implements the Wasserstein distance on (R,Delta_infinity).

**Syntax:**
```
%  res = uW_R(pos,mu1,mu2,p)
%   pos - joint support of mu1 and mu2
%   mu1 - probability vector that indicates the mass of mu1 on the corresponding position in the joint support
%   mu2 - probability vector that indicates the mass of mu2 on the corresponding position in the joint support
%   p   - real number >=1

%   
% Returns:
%   res   - the (R,Delta_infinity)-Wasserstein distance between the probability measures mu1 and mu2 to the power of p.

```
## uGW-distance
### ultraGWcgd.m 
**ultraGWcgd.m** is a Matlab function that approximates the ultrametric Gromov-Wasserstein distance of order p via conditional gradient descent (Frank-Wolfe algorithm with fixed decreasing step size). It is possible to start the gradient descent algorithm from multiple random couplings.

**Syntax:**
```
% [res,coupl,result_vec,result_coupl_cell] = ultraGWcgd(ux,uy,mux,muy,p,iterations,num_samples, num_skips)
%   ux  - ultrametric distance matrix
%   uy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%   iterations - number of iterations for the gradient descent
%   num_samples - number of random starting couplings for the gradient descent (if this is set to 0 only the independece coupling is used)
%   num_skips - number of skips between two random couplings 
%   
% Returns:
%   res   - the best approximation of the ultrametric Gromov-Wasserstein distance of order p between ux and uy with probability measures mux and muy, respectively.
%   plan  - the corresponding coupling
%   result_vec - a vector containing all stationary points of the gradient descent started from the different couplings
%   result_plan_cell - the collection of corresponding couplings
```


# dGW
For comparison, we also include codes for functions related to the usual Gromov-Wasserstein distance.
## Lower Bounds
### dgwslb.m 
**dgwslb.m** is a Matlab function for computing the second lower bound of dGW of two metric measure spaces.

**Syntax:**
```
% res = dgwslb(dx,dy,mux,muy,p)
%   dx  - ultrametric distance matrix
%   dy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res - pth-uGW-second lower bound between dx and dy with probability measures mux and muy, respectively.
```
### dgwtlb.m 
**dgwtlb.m** is a Matlab function for computing the third lower bound of dGW of two metric measure spaces.

**Syntax:**
```
% [res,coupl] = dgwtlb(dx,dy,mux,muy,p)
%   dx  - ultrametric distance matrix
%   dy  - ultrametric distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%
% Returns:
%   res    - pth-uGW-third lower bound between dx and dy with probability measures mux and muy, respectively.
%   coupl  - the corresponding optimal coupling
```
### dW_Rp.m
**dW_Rp.m** is a Matlab function that implements the Wasserstein distance of order p on the real line.

**Syntax:**
```
%  res = dW_Rp(pos1,pos2,mu1,mu2,p)
%   pos1  - support of mu1
%   pos2  - support of mu2
%   mu1 - probability vector
%   mu2 - probability vector
%   p   - real number >=1

%   
% Returns:
%   res   - the Wasserstein distance between the probability measures mu1 and mu2 to the power of p.

```
## dGW-distance
### dGWcgd.m 
**dGWcgd.m** is a Matlab function that approximates the Gromov-Wasserstein distance of order p via conditional gradient descent (Frank-Wolfe algorithm with fixed decreasing step size). It is possible to choose multtiple random starting points for the gradient descent algorithm.

**Syntax:**
```
% [res,coupl,result_vec,result_coupl_cell] = dGWcgd(dx,dy,mux,muy,p,iterations,num_samples, num_skips)
%   dx  - distance matrix
%   dy  - distance matrix
%   mux - probability vector
%   muy - probability vector
%   p   - real number >=1
%   iterations - number of iterations for the gradient descent
%   num_samples - number of random starting couplings for the gradient descent (if this is set to 0 only the independece coupling is used)
%   num_skips - number of skips between two random couplings 
%   
% Returns:
%   res   - the best approximation of the Gromov-Wasserstein distance of order p between dx and dy with probability measures mux and muy, respectively.
%   plan  - the corresponding coupling
%   result_vec - a vector containing all stationary points of the gradient descent started from the different couplings
%   result_plan_cell - the collection of corresponding couplings
```
## Tests
In this file, we includes codes for computing uGW-slb and dGW-slb between randomly generated ultrametric spaces.
### uGW-test.m 
**uGW-test.m** is a Matlabscript for computing uSLB,uTLB and uGW between each pair of ultrametric spaces from a given file **Us.mat** with uniform probability measures and generates corresponding distance matrices sorted by number of blocks in ultrametric spaces. Before running the code, one should first run either one of the following Matlab codes to generate random ultrametric spaces recorded in file **Us.mat**:

* **test_subsample.m** Subsample ultrametric spaces from some given randomly generated ultrametric spaces.
* **test_independent.m** Generates ultrametric spaces independently.
* 
### dGW-test.m 
**dGW-test.m** is a Matlabscript for computing dSLB,dTLB and dGW between each pair of ultrametric spaces from a given file **Us.mat** with uniform probability measures and generates corresponding distance matrices sorted by number of blocks in ultrametric spaces. Before running the code, one should first run either one of the following Matlab codes to generate random ultrametric spaces recorded in file **Us.mat**:

* **test_subsample.m** Subsample ultrametric spaces from some given randomly generated ultrametric spaces.
* **test_independent.m** Generates ultrametric spaces independently.

### visual.m 
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
