#generate data
rm(list=ls())
library(MASS)
N = 1000 #Sample size
mu = rep(0,3)
rho = 0.5

Sigma = matrix(c(1,rho,rho^2,rho,1,rho,rho^2,rho,1),ncol = 3)
X = mvrnorm(N,mu, Sigma)# Covariates

lambdaC = 0.5 #Censoring time C~exponential(lambdaC) 
C = rexp(N,lambdaC)



# beta0_true = runif(1,-3,3)
beta0_true = 0
beta_true = matrix(runif(length(mu),-3,3),ncol = 1)
error = rnorm(N,0,1)
g1 = function(u,gamma0){
  return(exp(-gamma0 * u))
}

# log_T = beta0_true + g1(X %*% beta_true, gamma0 = 1.4) + error
# log_T = (X %*% beta_true) + error
# #
# log_Y = pmin(log_T,log(C))
# Delta = log_T - log_Y
# Delta[which(Delta>0)]=1

log_Y = (X %*% beta_true) + error
Delta = rbinom(N,1,0.6)
