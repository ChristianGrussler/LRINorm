# LRINorm
MATLAB code for Low-rank Optimization by Low-Rank Inducing Norms and their non-convex counterparts. 

## Publications:

* [Rank Reduction with Convex Constraints](https://lup.lub.lu.se/search/publication/54cb814f-59fe-4bc9-a7ef-773cbcf06889)
* [Low-rank Inducing Norms with Optimality Interpretations](https://arxiv.org/abs/1612.03186)
* [Low-rank Optimization with Convex Constraints](https://arxiv.org/abs/1606.01793)
* [Local Convergence of Proximal Splittinge Methods for Rank Constrained Problems](https://arxiv.org/abs/1710.04248)
* [The Use of the r* Heuristic in Covariance Completion Problems](http://www.control.lth.se/index.php?mact=ReglerPublicationsB,cntnt01,showpublication,0&cntnt01LUPid=a61669c7-29b9-41ee-82da-9c825b08f8d8&cntnt01returnid=60)
* [On optimal low-rank approximation of non-negative matirces](http://lup.lub.lu.se/search/ws/files/21812505/2015cdcGrusslerRantzer.pdf)

## Prerequisites
No prerequisites for proximal splitting methods. 
For CVX implementations get CVX from: http://cvxr.com/cvx/

## Installation

Download all files and add to path: https://github.com/LowRankOpt/LRINorm/archive/master.zip 

## Documentation

### Examples:
There are three examples in the "Example" folder:

1. Exact Matrix Completion
2. Low-rank approximation with Hankel constraint
3. Model Order Reduction through Hankel low-rank approximation


### Proximal Mappings:
The folder "Prox" contains the proximal mappings to the low-rank inducing Frobenius and Spectral norm as well as their non-convex counter parts.

For low-rank inducing Frobenius norm: p = 2
For low-rank inducing Frobenius norm: p = 1

#### Low-rank inducing Spectral and Frobenius norms: 

Proximal mapping of the low-rank inducing norms with at Z with parameter r and scaling factor gamma:
```
X = proxnormrast(Z,r,p,gamma)
```
#### Squared Low-rank inducing Spectral and Frobenius norms: 
Proximal mapping of the SQUARED low-rank inducing norms with at Z with parameter r and scaling factor gamma :
```
X = proxnormrast_square(Z,r,p,gamma)
```
#### Projection on the epi-graph of the low-rank inducing norms: 
Projection of (Z,zv) on the epi-graph of the low-rank inducing norms with parameter r and scaling factor gamma:
```
[X,xv] = proxnormrast_square(Z,zv,r,p,gamma)
```

#### Non-convex proximal mappings for Frobenius and Spectral norm: 

Non-convex proximal mapping of at Z with parameter r and scaling factor gamma:
```
X = proxnonconv(Z,r,p,gamma)
```
#### Non-convex proximal mappings for squared Frobenius and Spectral norm:
Non-convex proximal mapping for the SQUARED norms at Z with parameter r and scaling factor gamma :
```
X = proxnonconv_square(Z,r,p,gamma)


