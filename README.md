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

## ugwtlb.m 
**ugwtlb.m** is a Matlab function computing the third lower bound of uGW of two ultrametric measure spaces.

**Syntax:**
```
% d = ugwtlb(u1,u2,mu1,mu2,p)
%   u1  - ultrametric distance matrix
%   u2  - ultrametric distance matrix
%   mu1 - probability vector
%   mu2 - probability vector
%   p   - real number >=1
%
% Returns:
%   res   - pth-uGW-third lower bound between u1 and u2 with probability measures mu1 and mu2, respectively.
%   plan  - the corresponding transport plan
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
