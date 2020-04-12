# uGW

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
