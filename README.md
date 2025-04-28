# R package 'rlpls'

This is a package for sparse partial least squares regression with random LASSO (RLPLS). RLPLS is a statistical framework derived from sparse partial least squares (SPLS) and random LASSO methodologies.

# Installation
To install the package, it requires the package 'devtools'.
```
install.packages("devtools")
devtools::install_github("jinky00/rlpls")
```

# Example
```
set.seed(1111)
# Set the size of samples and true coefficients
n = 30
beta0 = c(runif(20,-10,-5),runif(20,5,10),rep(0,10))

# Construct predictors & response
# the structure is similar with that from Chun and Keles (2010)
xxx = yyy = c(NULL)
H1 = c(rep(2,times=25),rep(4,times=n-25))
H2 = rep(10,n)
for(i in 1:length(beta0)){
  if(i <= 30){
    xxx = cbind(xxx, H1+rnorm(n,0,1))
  }else{
    xxx = cbind(xxx, H2+rnorm(n,0,1))
  }
};rm(i)# for i

yyy = xxx %*% beta0 + rnorm(n,0,1.5^2)

# Run RLPLS
library(rlpls); library(dplyr)
result = rlpls(X = xxx, Y = yyy, b_var_num = c(25,20), Btimes = 5000 , variable_select = "simpls")
```

# References
In Preparation
