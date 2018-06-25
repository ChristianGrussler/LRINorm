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

Serveral examples are found in the folder "Examples", which demonstrate how to use the other files.

The folder "Prox" contains several proximal mappings to the low-rank inducing Frobenius and Spectral norm:

### Low-rank inducing Frobenius norm: 

Proximal mapping of the low-rank inducing Frobenius norm with parameter r and scaling factor gamma is evaluated at a matrix Z by:
```
X = proxnormrast(Z,r,2,gamma)
```
Proximal mappings of the SQUARED low-rank inducing Frobenius norm with parameter r and scaling factor gamma is evaluated at a matrix Z by:
```
X = proxnormrast_square(Z,r,2,gamma)
```
Projection on the ep-graph of the low-rank inducing Frobenius norm with parameter r and scaling factor gamma is evaluated at (Z,zv) by:
```
[X,xv] = proxnormrast_square(Z,zv,r,2,gamma)
```
### Low-rank inducing Spectral norm: 

Proximal mapping of the low-rank inducing Frobenius norm with parameter r and scaling factor gamma is evaluated at a matrix Z by:
```
X = proxnormrast(Z,r,1,gamma)
```
Proximal mappings of the SQUARED low-rank inducing Frobenius norm with parameter r and scaling factor gamma is evaluated at a matrix Z by:
```
X = proxnormrast_square(Z,r,1,gamma)
```
Projection on the ep-graph of the low-rank inducing Frobenius norm with parameter r and scaling factor gamma is evaluated at (Z,zv) by:
```
[X,xv] = proxnormrast_square(Z,zv,r,1,gamma)
```

