# LRINorm
MATLAB code for Low-rank Optimization by Low-Rank Inducing Norms and non-convex Douglas-Rachford. 

## Purpose:
Low-rank rank inducing norms and non-convex Proximal Splitting Algoriths attempt to find exact rank/cardinality-r solutions to minimization problems with convex loss functions, i.e., avoiding of regularzation heuristics. This project provides MATLAB implementations for the proximal mappings of the low-rank inducing Frobenius and Spectral norms, as well as, their epi-graph projections and non-convex counter parts.

## Literature:

### Low-rank inducing norms: 
* [Rank Reduction with Convex Constraints](https://lup.lub.lu.se/search/publication/54cb814f-59fe-4bc9-a7ef-773cbcf06889)
* [Low-rank Inducing Norms with Optimality Interpretations](https://arxiv.org/abs/1612.03186)
* [Low-rank Optimization with Convex Constraints](https://arxiv.org/abs/1606.01793)
* [The Use of the r* Heuristic in Covariance Completion Problems](http://www.control.lth.se/index.php?mact=ReglerPublicationsB,cntnt01,showpublication,0&cntnt01LUPid=a61669c7-29b9-41ee-82da-9c825b08f8d8&cntnt01returnid=60)
* [On optimal low-rank approximation of non-negative matirces](http://lup.lub.lu.se/search/ws/files/21812505/2015cdcGrusslerRantzer.pdf)

### Non-convex counter parts:
* [Local Convergence of Proximal Splittinge Methods for Rank Constrained Problems](https://arxiv.org/abs/1710.04248)

## Prerequisites
No prerequisites for proximal splitting methods. 
For CVX implementations get CVX from: http://cvxr.com/cvx/

## Installation

Download all files and add to path: https://github.com/LowRankOpt/LRINorm/archive/master.zip 

## Documentation
In the following it holds that
* for the low-rank inducing Frobenius norm: p = 2
* for the low-rank inducing Spectral norm:  p = inf

### Examples
There are three examples in the "Example" folder:

1. Exact Matrix Completion
2. Low-rank approximation with Hankel constraint
3. Model Order Reduction through Hankel low-rank approximation

### Optimization

The "Optimizaton" folder contains Douglas-Rachford splitting as well as CVX implementations for the low-rank inducing Frobenius and Spectral norms. It is easy to modify these functions for other constraints!

#### Exact Matrix completion

Let N be a matrix and Index be a binary matrix of the same size, where the ones indicate the known entries N. We attempt to find a rank-r completion M through

i) Low-rank inducing norms with Douglas-Rachford splitting:

```
M = drcomplete(N,Index,r,p)
```

ii) Low-rank inducing norms with CVX:

```
M = cvxcomplete(N,Index,r,p)
```

iii) Non-convex Douglas-Rachford splitting:

```
M = drcomplete(N,Index,r,p,'solver','NDR')
```

#### Low-rank Hankel Approximation

Let H be a matrix. We attempt to find a rank-r Hankel approximation M that minimizes the Frobenius norm error through

i) Low-rank inducing Frobnius norms with Douglas-Rachford splitting:

```
M = drhankelapprox(H,r)
```

ii) Low-rank inducing Frobenius norms with CVX:

```
M = cvxhankelapprox(H,r)
```

iii) Non-convex Douglas-Rachford splitting with Frobenius norm:

```
M = drhankelapprox(H,r,'solver','NDR')
```

### Proximal Mappings
The folder "Prox" contains the proximal mappings to the low-rank inducing Frobenius and Spectral norm, as well as, as their epi-graph projections and non-convex counter parts.

#### Low-rank inducing Spectral and Frobenius norms: 

Proximal mapping of the low-rank inducing norms with at Z with parameter r and scaling factor gamma:
```
X = proxnormrast(Z,r,p,gamma)
```
#### Squared Low-rank inducing Spectral and Frobenius norms: 
Proximal mapping of the SQUARED low-rank inducing norms with at Z with parameter r and scaling factor gamma:
```
X = proxnormrast_square(Z,r,p,gamma)
```
#### Projection on the epi-graph of the low-rank inducing norms: 
Projection of (Z,zv) on the epi-graph of the low-rank inducing norms with parameter r and scaling factor gamma:
```
[X,xv] = projrast(Z,zv,r,p,gamma)
```

#### Non-convex proximal mappings for Frobenius and Spectral norm: 

Non-convex proximal mapping of at Z with parameter r and scaling factor gamma:
```
X = proxnonconv(Z,r,p,gamma)
```
#### Non-convex proximal mappings for squared Frobenius and Spectral norm:
Non-convex proximal mapping for the SQUARED norms at Z with parameter r and scaling factor gamma:
```
X = proxnonconv_square(Z,r,p,gamma)
```
