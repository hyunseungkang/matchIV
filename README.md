# matchIV
matchIV contains a set of R functions to do matching-based IV analysis with a single binary IV. The details of these methods are in Kang, Kreuels, May, and Small (2016) and Kang, Peck and Keele (2016).

## Installation

To load these functions into R, run the following commands
```R
source("https://raw.githubusercontent.com/hyunseungkang/invalidIV/master/TSHT.R")
```

## Examples
The code example.r has additional working examples.

```R
### Obtain the two-stage hard thresholding (TSHT) code from Github ###	
source("https://raw.githubusercontent.com/hyunseungkang/matchIV/master/effectRatio.r")

### R Packages to Load ###
# The AER and MASS packages are only needed to run the working example.
library(AER)
library(MASS)

### Working Example ###
### n = 500, pz = 10 IVs (s = 3 invalid)
# Y: n by 1 vector of outcomes (must be continuous)
# D: n by 1 vector of treatments (continuous or discrete)
# Z: n by pz vector of instruments (continuous or discrete)
# beta:	  true treatment effect (set at 1)

# Create data #
library(MASS)
n = 500; L = 10; s = 3
alpha = c(rep(3,s),rep(0,L-s)); beta = 1; gamma = c(rep(1,L))
epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
Z = matrix(rnorm(n*L),n,L)

epsilon = mvrnorm(n,rep(0,2),epsilonSigma)
D = Z %*% gamma + epsilon[,1]
Y = Z %*% alpha + D * beta + epsilon[,2]

### Oracle Two-Stage Least Squares Estimator ###
# This oracle knows exactly which IVs are invalid
# and should perform very well,	  with the point estimate
# associated with D close to 1 and the 95% confidence
# interval covering 1

summary(ivreg(Y ~ D + Z[,1:s] - 1 | Z - 1))
confint(ivreg(Y ~ D + Z[,1:s] - 1 | Z - 1))[1,]


### Our TSHT Estimator ###
# Our estimator	does not assume	which IVs are invalid
# a priori. It should perform as well as the the oracle
# above	as sample size increases.

# Output is a list that includes
# beta: point estimate of the treatment effect
# se: standard error of the beta
# ci: 1 - alpha confidence interval for beta
# V: estimated set of valid IVs
# S: estimated set of relevant IVs
 
TSHT.ldim(Y,D,Z)
```

## References 
Guo, Z., Kang, H., Cai, T. T., and Small, D. S. (2016). <a href="http://arxiv.org/abs/1603.05224">Confidence Interval for Causal Effects with Invalid Instruments using Two-Stage Hard Thresholding.</a> Technical Report.

Kang, H., Cai, T. T., Small, D. S. (2016). <a href="http://arxiv.org/abs/1504.03718">A Simple and Robust Confidence Interval for Causal Effects with Possibly Invalid Instruments.</a> Technical Report.
