### TITLE: Full Matching Effect Ratio Code ###
### MAINTAINER: Hyunseung Kang (hskang at stanford dot edu)
### LAST UPDATE: Nov. 29, 2015
### SUMMARY: The code implements the full matching IV method 
###          described in the paper Kang, Kreuels, May, 
###          and Small (2015+). Specifically, the main 
###          function, estLambda(), estimates the 
###          effect ratio parameter described in the 
###          said paper.
### REQUIRES: MASS, optmatch, AER (R packages, install them 
###           from CRAN) 
### USAGE: Copy-paste the R code below. At the end of
###        the code is a simple running example of
###        the code from simulated data. 

### Libraries to Load ####
library(optmatch)
library(AER)

### Functions ###
##### smahal: Obtain rank-based Mahalanobis distance
### FUNCTION: Creates a rank-based Mahalanobis distance. 
### INPUT: 1) Z: instruments (n by 1 vector), Z =1 (treated) and Z = 0 (not treated)
###        2) X: covariates  (n by p matrix}
### OUTPUT: 1) n by n rank-based Mahalanobis distance matrix 
smahal= function(Z,X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  rownames(X)<-1:n
  k<-dim(X)[2]
  m<-sum(Z)
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,m,n-m)
  Xc<-X[Z==0,]
  Xt<-X[Z==1,]
  rownames(out)<-rownames(X)[Z==1]
  colnames(out)<-rownames(X)[Z==0]
  icov<-ginv(cv)
  for (i in 1:m) out[i,]<-mahalanobis(Xc,Xt[i,],icov,inverted=T)
  out
}

##### addcaliper: adds caliper to the matching matrix
### FUNCTION: Creates a calipered rank-based Mahalanobis distance from a distance matrix.
### INPUT: 1) dmat: distance matrix.
###        2) Z: instruments (n by 1 vector), Z =1 (treated) and Z = 0 (not treated)
###        3) logitp: logit propensity scores from logistic regression (n by 1 vector)
###        4) calipersd, indicator on when to caliper
###        5) penalty, how much penalty for exceeding caliper bounds
### OUTPUT: 1) n by n with caliper added
addcaliper=function(dmat,Z,logitp,calipersd=.2,penalty=1000){
  sd.logitp=sd(logitp)
  adif=abs(outer(logitp[Z==1],logitp[Z==0],"-"))
  adif=(adif-(calipersd*sd.logitp))*(adif>(calipersd*sd.logitp))
  dmat=dmat+adif*penalty
  dmat
}

##### matching: does full matching
### FUNCTION: creates a vector that matches treatment to control
### INPUT: 1) Z: a n-length vector that lists treatment (as 1) and control (as 0)
###        2) X: a n-by-p matrix of covariates
###        3) type: type of matching, defaults to full matching
###        4) ratio: if full matching or 1:m matching, checks to see the control/treatment ratio
###        5) matchInfo: should the function print out the matching information?
### OUTPUT: 1) a n-vector vector containing match numbers.
### NOTE: The matching function (smahal, prop.score, caliper, and fullmatch) are
###       very sensitive to how the observations are ordered. Currentyl, smahal
###       orders the observations by how Z is ordered.
matching = function(Z,X,type=c("full","pair","1:m"),ratio,matchInfo=FALSE) {
  type = match.arg(type)
  
  # Pre-processing
  Z = as.logical(Z);
  
  # Rank-based distance matrix
  rank.mat = smahal(Z,X)
  
  # Propensity scores
  model = glm(Z ~ X,family="binomial"); prop.score = predict(model);
  
  # Fix rank-based distance matrix by propensity score calipers
  rank.mat.caliper = addcaliper(rank.mat,Z,prop.score)

  # Actual matching
  if(type == "pair") pairmatchvec = pairmatch(rank.mat.caliper,data=1:length(Z))
  if(type == "full") {
  	  if(missing(ratio)) pairmatchvec = fullmatch(rank.mat.caliper,data=1:length(Z))
  	  else pairmatchvec = fullmatch(rank.mat.caliper,data=1:length(Z),min.controls=1/ratio,max.controls=ratio)
  }
  if(type == "1:m") pairmatchvec = pairmatch(rank.mat.caliper,data=1:length(Z),controls=ratio)
  if(matchInfo) print(summary(pairmatchvec))
  
  # Post-processing
  pairmatchvec = as.numeric(substr(pairmatchvec,start=3,stop=20))  
  return(pairmatchvec) 
}

##### estLambda: estimates Lambda, provide confidence interval, and p-value
### FUNCTION: provides point est, confidence interval, and p-value
### INPUT: Z: instrument (n by 1 vector)
###        R: response (n by 1 vector)
###        D: dose (n by 1 vector)
###        X: covariate matrix (n by p) 
###        matchedNumber: matching vector from matching() function
###        null: the lambda_0 value in H0: lambda = lambda_0
###        alphaLevel: alpha level for the confidence interval
### OUTPUT: point Est, confidence interval, and p-value
estLambda = function(Z,R,D,X,matchedNumber,null = 0,alphaLevel = 0.05) {
  if(alphaLevel >= 0.5) stop("alphaLevel is improperly set")
  
  X.factorize = model.matrix(~X) 
  if(missing(matchedNumber)) matchedNumber = matching(Z,X,"full")
  			     
  # Deal with non-matched individuals #
  R.matchedIndiv = R[!is.na(matchedNumber)]
  D.matchedIndiv = D[!is.na(matchedNumber)]
  Z.matchedIndiv = Z[!is.na(matchedNumber)]
  matchedNumber.matchedIndiv = matchedNumber[!is.na(matchedNumber)]

  # Sort individuals in ascending matching order
  matchedNumber.sortedIndex = sort(matchedNumber.matchedIndiv,index.return=TRUE)$ix
  R.sorted = R.matchedIndiv[matchedNumber.sortedIndex]
  D.sorted = D.matchedIndiv[matchedNumber.sortedIndex]
  Z.sorted = Z.matchedIndiv[matchedNumber.sortedIndex]
  matchedNumber.sorted = matchedNumber.matchedIndiv[matchedNumber.sortedIndex]
  
  # Calculate the size of each matched set and the corresponding weights
  ni = tabulate(matchedNumber.sorted); 
  ni = ni[!(ni == 0)] #this cleans up a flaw with matching algorithm 
                      #where some matched numbers (e.g. 1:I) are not used between
		      	     	  	  	    #1 to I
  I = length(ni)
  wi = ni^2 / (ni - 1)  #Actual formula: ni^2 /(mi * (ni - mi))

  # Calculate Vi, Gi, and Hi
  # Formula
  # Vi = Gi - lambda * Hi
  # Gi = wi*(Zij - Zimean)(Rij - Rimean) = wi*(ZijRij - ni*Zmean.*Rmean.)
  # Hi = wi*(Zij - Zimean)(Dij - Dimean) = wi*(ZijDij - ni*Zmean.*Dmean.)

  # This requires some thought. Basically, because R,D,Z are all ordered based on matchedset number,
  # we use cumsum to obtain means of R, D, and Z within each group as well as means of RZ and DZ.
  # The key here is using the cumsum(ni), which keeps track of the indices in all the ordered vectors
  ni.cumsum = cumsum(ni)
  RZ.cumsum = cumsum(R.sorted*Z.sorted); DZ.cumsum = cumsum(D.sorted*Z.sorted)
  R.cumsum = cumsum(R.sorted); D.cumsum = cumsum(D.sorted); Z.cumsum = cumsum(Z.sorted)
  R.means = (R.cumsum[ni.cumsum] - c(0,R.cumsum[ni.cumsum[-length(ni.cumsum)]]))/ni
  Z.means = (Z.cumsum[ni.cumsum] - c(0,Z.cumsum[ni.cumsum[-length(ni.cumsum)]]))/ni
  D.means = (D.cumsum[ni.cumsum] - c(0,D.cumsum[ni.cumsum[-length(ni.cumsum)]]))/ni
  
  RZ.group = RZ.cumsum[ni.cumsum] - c(0,RZ.cumsum[ni.cumsum[-length(ni.cumsum)]])
  DZ.group = DZ.cumsum[ni.cumsum] - c(0,DZ.cumsum[ni.cumsum[-length(ni.cumsum)]])

  Gi = wi * RZ.group - wi * ni * R.means * Z.means
  Hi = wi * DZ.group - wi * ni * D.means * Z.means
  
  # Compute the point estimate 
  pointEst = sum(Gi) / sum(Hi)
  
  # Compute the p-value under the null
  ViNull = Gi - null * Hi
  testStatNull = mean(ViNull) / sqrt(1/(I *(I-1)) * sum( (ViNull - mean(ViNull))^2) )
  pvalue = 2* (1 - pnorm(abs(testStatNull)))
  
  # Compute the quadratic terms for CI
  # A*lambda^2 + B*lambda + C = 0 (technically, it's A*lambda^2 + B*lambda + C <= 0)
  q = qnorm(1 - alphaLevel/2)
  A = sum(Hi)^2/I^2 - q^2 /(I * (I-1)) * sum( (Hi - mean(Hi))^2)
  B = -2 * (sum(Hi) * sum(Gi) / I^2 - q^2 / (I * (I-1)) * sum( (Hi - mean(Hi)) * (Gi - mean(Gi))))
  C = sum(Gi)^2/I^2 - q^2/(I * (I-1)) * sum( (Gi - mean(Gi))^2)

  detQuad = round(B^2 - 4*A*C,9) #6 is set for numerical accuracy
  if( detQuad <= 0) {
    if(A < 0) {
      cis = matrix(c(-Inf,Inf),1,2)
      } else {
        cis = matrix(c(NA,NA),1,2)
	} 
  }
  if(detQuad > 0) {
    if(A < 0) {
      up.ci = (-B - sqrt(detQuad))/ (2*A) 
      low.ci = (-B + sqrt(detQuad)) / (2*A) 
        cis = matrix(c(-Inf,low.ci,up.ci,Inf),2,2,byrow=TRUE)
	} else {
	  low.ci = (-B - sqrt(detQuad))/ (2*A) 
      up.ci = (-B + sqrt(detQuad)) / (2*A) 
        cis = matrix(c(low.ci,up.ci),1,2)
	}
  }
  return(list(pointEst = pointEst,cis = cis,pvalue = pvalue))
} 
###### est2SLS
### FUNCTION: computes the 2SLS estimate using the AER package
### INPUT: Z: instrument (n by 1 vector)
###        R: response (n by 1 vector)
###        D: dose (n by 1 vector)
###        X: covariate matrix (n by p) 
###        null: null Value under H_0:
### OUTPUT: point Est, confidence interval, and p-value
est2SLS <- function(Z,R,D,X,null=0) {
  library(AER)
  model = ivreg(R ~ D + X | Z + X)
  model.summary = summary(model)
  coefficientTable = model.summary$coefficients
  dfTest = model.summary$df[2]
  estValue = as.numeric(coefficientTable[2,1])
  stdError = as.numeric(coefficientTable[2,2])
  tValue = abs( (estValue - null) / stdError)
  pvalue = 2*(1 - pt(tValue,dfTest))
  confInt = matrix(as.numeric(confint(model)[2,]),1,2)
  return(list(pointEst = estValue,cis = confInt,pvalue = pvalue))
}

###### estOLS
### FUNCTION: computes the OLS estimate 
### INPUT: Z: instrument (n by 1 vector), note this is not used in OLS
###        R: response (n by 1 vector)
###        D: dose (n by 1 vector)
###        X: covariate matrix (n by p) 
### OUTPUT: point Est, confidence interval, and p-value
estOLS <- function(Z,R,D,X) {
  model = lm(R ~ D + X)
  estValue = as.numeric(coefficients(model))[2]
  pvalue = as.numeric(summary(model)$coefficients[2,4])
  confInt = matrix(as.numeric(confint(model)[2,]),1,2)
  return(list(pointEst = estValue,cis = confInt,pvalue=pvalue))
}

##### sensInt 
### FUNCTION: does sensitivity analysis
### INPUT: Gamma: a scalar value 
###        Z: binary treatment/control vector (n time 1) 
###        R: binary response (n * 1)
###        pairmatchvec: vector of matches (handles no matches)
### OUTPUT: lower and upper bounds for sensitivity analysis
### NOTE: Can only handle sensitivity analysis if there are multiple controls 
###      (but not if there are multiple treatments) AND if the response is binary
sensInt <- function(Gamma,Z,R,pairmatchvec) {
  matchedSets = unique(pairmatchvec[!is.na(pairmatchvec)]); 
  cs.plus = rep(0,length(I))
  ps.plus = rep(0,length(I))
  ps.minus = rep(0,length(I))
  testStat = rep(0,length(I))
  ds= rep(0,length(I))
  for(i in 1:length(matchedSets)) { 
    index = matchedSets[i]
    matchedSetID = which(pairmatchvec == index)
    ni = length(matchedSetID); mi = sum(Z[matchedSetID])	  
    cs.plus[i] = sum(R[matchedSetID])
    ps.plus[i] = Gamma * cs.plus[i]/(Gamma*cs.plus[i] + ni - cs.plus[i])
    ps.minus[i] = cs.plus[i]/(cs.plus[i] + (ni - cs.plus[i])*Gamma)
    ds[i] = ni/(mi*(ni - mi))
    testStat[i] = ds[i] * sum(Z[matchedSetID] * R[matchedSetID])
  }
  T = sum(testStat)
  lowerBound = (T - sum(ds * ps.plus)) /(sqrt(sum(ds^2 * ps.plus * (1 - ps.plus))))
  upperBound = (T - sum(ds * ps.minus)) /(sqrt(sum(ds^2 * ps.minus * (1 - ps.minus))))
  return(c(lowerBound,upperBound))
}

##### balanceCheck 
### FUNCTION: checks the balance of covariates before and after matching
### INPUT: Z: treatment/control vector (n time 1) 
###        X: covariates (n * p)
###        pairmatchvec: vector of matches (handles no matches)
### OUTPUT: a table of standardized differences
balanceCheck <- function(Z,X,pairmatchvec) {
   # Before Matching #
   X.C.BeforeMatch = X[Z == 0,]; X.T.BeforeMatch = X[Z == 1,];
   overallSD = sqrt(apply(X.C.BeforeMatch,2,var) + apply(X.T.BeforeMatch,2,var)/2)
   stdDiff.Before = abs((colMeans(X.T.BeforeMatch) - colMeans(X.C.BeforeMatch))/overallSD)
   # After Matching #
   matchedSets = unique(pairmatchvec[!is.na(pairmatchvec)]); 
   X.C.AfterMatch = matrix(0,length(matchedSets),ncol(X)); X.T.AfterMatch = matrix(0,length(matchedSets),ncol(X));
   nTreatedAfterMatch = rep(0,length(matchedSets))
   for(i in 1:length(matchedSets)) { 
      treatedZ = which(pairmatchvec == matchedSets[i] & Z == 1)
     controlZ = which(pairmatchvec == matchedSets[i] & Z == 0)
     X.C.AfterMatch[i,] = colMeans(X[controlZ,,drop=FALSE]); X.T.AfterMatch[i,] = colMeans(X[treatedZ,,drop=FALSE])
     nTreatedAfterMatch[i] = length(treatedZ)
   }
   stdDiff.After = apply( X.T.AfterMatch - X.C.AfterMatch,2,function(x){ sum(x * nTreatedAfterMatch) / sum(nTreatedAfterMatch)}) / overallSD
   stdDiff.After = abs(stdDiff.After)
   output = cbind(stdDiff.Before,stdDiff.After)
   colnames(output) = c("Standardized differences before matching","Standardized differences after matching")
   return(output)
}

##### effectSampleSize 
### FUNCTION: computes the asymptotic variance of effectRatio estimator and the maximum strata size for effectratio estimator
### INPUT: Z: treatment/control vector (n time 1) 
###        matchedNumber: vector of matches
### OUTPUT: a two-dimensional vector, the first dimension being the effective sample size and the second dimension being the maximum strata size.
effectSampleSize <- function(Z,matchedNumber) {
  # Deal with non-matched individuals #
  Z.matchedIndiv = Z[!is.na(matchedNumber)]
  matchedNumber.matchedIndiv = matchedNumber[!is.na(matchedNumber)]

  # Sort individuals in ascending matching order
  matchedNumber.sortedIndex = sort(matchedNumber.matchedIndiv,index.return=TRUE)$ix
  Z.sorted = Z.matchedIndiv[matchedNumber.sortedIndex]
  matchedNumber.sorted = matchedNumber.matchedIndiv[matchedNumber.sortedIndex]
  
  # Calculate the size of each matched set and the corresponding weights
  ni = tabulate(matchedNumber.sorted); I = length(ni)
  ni.max = max(ni)
  asympVarEffectRatio = (sum(ni^3 / (ni - 1)) * 1/I) / (sum(ni) / I)^2
  print(paste("Asymptotic variance of effect ratio estimator is",asympVarEffectRatio))
  print(paste("Maximmum strata size is",ni.max))
  return(c(asympVarEffectRatio,ni.max))
}

##### matchSizeAnalysis
### FUNCTION: analyzes the size of matched Sets as a function of the covariates
### INPUT: pairmatchvec: vector of of matches
###        Z: treatment/control vector 
###        X: covariates
### OUTPUT: a two-dimensional vector, the first dimension being the effective sample size and the second dimension being the maximum strata size.
matchSizeAnalysis <- function(pairmatchvec,Z,X) {
  matchedSets = unique(pairmatchvec[!is.na(pairmatchvec)]); I = length(matchedSets)
  Y.sizeMatch = rep(0,I)
  Y.controlMatch = rep(0,I)
  X.covAvg = matrix(0,I,ncol(X))
  
  for(i in 1:I) {
    set.i = which(matchedSets[i] == pairmatchvec)
    Y.controlMatch[i] = sum(Z[set.i])
    Y.sizeMatch[i] = length(set.i)
    for(j in 1:ncol(X.covAvg)) {
      X.covAvg[i,j] = mean(X[set.i,j])
      }
  }
  colnames(X.covAvg) = colnames(X)
  return(data.frame(Y.sizeMatch = Y.sizeMatch,Y.controlMatch=Y.controlMatch,X.covAvg))
}
