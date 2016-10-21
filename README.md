# matchIV
matchIV contains a set of R functions to do matching-based IV analysis with a single binary IV. The details of these methods are in Kang, Kreuels, May, and Small (2016) and Kang, Peck and Keele (2016).

## Installation

To load these functions into R, run the following commands
```R
source("https://raw.githubusercontent.com/hyunseungkang/matchIV/master/effectRatio.r")
```

## Examples
The code example.r has additional working examples.

```R
### Obtain effectRatio code from Github ###	
source("https://raw.githubusercontent.com/hyunseungkang/matchIV/master/effectRatio.r")

##### Example Run #####
### Generate simulated data ###

### Simulated data contains n=1000 individuals, 
### the response (denoted as R), the exposure
### (denoted as D), the instrument (denoted as
###  Z), and the covariates (denoted as X)

### Simulated Data Characteristics
### R is continuous (e.g. height)
### D is continuous (e.g. amount of malaria parasite in plasma)
### Z is binary (1/0 e.g. Sickle cell trait status)
### X1 is a continuous variable (e.g. birthweight),
### X2 is a binary variable (1/0, e.g. sex)
### NOTE1: our method can handle continuous, discrete, binary R,D,Xs, 
###       and many covariates. Our method only
###       requires that Z is binary (1/0) and non-missing.
### NOTE2: If some values of Xs are missing 
###        (e.g. contains NAs), then
###        you must provide a "cleaned" X
###        without NAs, depending on how 
###        you want to take care of missing
###        values (e.g. imputation, MCAR, MAR, etc.)    
### NOTE3: If R or D are missing, then you must provide a 
###        cleaned R and D without NAs.

set.seed(1) #Random number generator seed. Change it as you please!
n = 1000
Z = runif(n) < 0.5
X1 = rnorm(n,rep(c(0,1),n/2))
X2 = rep(c(0,1),n/2)
X = cbind(X1,X2)

# Generate R and D
# Note that R and D are non-linear functions of X.
# -3 is the true value of the causal effect of D on R in this example
# Our point estimator should be close to -3 and our confidence interval
# should have the desired level of coverage
D = -1 + 0.2 * Z + X1 - 2*X1^2 - 2 * X2 + rnorm(n)
R = 2 - 3 * D + X1 - 2*X1^2 - 2 * X2 + rnorm(n)


### Run our method (automatically matches and estimateds Lambda)
estLambda(Z,R,D,X)

### Run 2SLS
est2SLS(Z,R,D,X)

### EXPECTED RESULT: 
###         In most cases, you will notice that our method (estLambda)'s
###         point estimate will be closer to -3. Furthermore, our confidence
###         interval should have the desired 95% level of coverage. 
###         In contrast, 2SLS (est2SLS) will be off from -3 most of the 
###         time and its confidence interval will not have the desired 
###         95% level of coverage.


#### Some other functions
### Match analysis
fullmatchvec = matching(Z,X) #full-match individuals
balanceCheck(Z,X,fullmatchvec) #checks balance of covariates
```

## References 
Kang, H., Peck, L., Keele, L. (2016). <a href="http://arxiv.org/abs/1606.04146">A Comparison of Inferential Techniques for Instrumental Variables Methods.</a> Technical Report.

Kang, H. (2016). <a href="http://journals.lww.com/epidem/Citation/publishahead/Matched_Instrumental_Variables___A_Possible.99005.aspx">Commentary: Matched Instrumental Variables: A Possible Solution to Severe Confounding in Matched Observational Studies?</a> <i> Epidemiology</i>,27, 624-632.

Kang, H., Kreuels, B., May, J., Small, D. S. (2016). <a href="https://projecteuclid.org/euclid.aoas/1458909919">Full Matching Approach to Instrumental Variables Estimation with Application to the Effect of Malaria on Stunting.</a> <i> Annals of Applied Statistics</i>,10,335-364.

