setwd("D:/2018+/term2/SIM/mycode")
source('GenerateData.r')
library(np)
# No Zero Denominator, used in C code for kernel estimation... (from 'np' package : https://CRAN.R-project.org/package=np)
NZD <- function(a) {
  sapply(1:NROW(a), function(i) {if(a[i] < 0) min(-.Machine$double.eps,a[i]) else max(.Machine$double.eps,a[i])})
}


# beta = matrix(runif(length(mu)+1,-3,3),ncol = 1)
# Xquta = cbind(matrix(rep(1,N),ncol = 1),X)

# beta0 = 0 no intercept
beta = matrix(runif(length(mu),-3,3),ncol = 1)
beta = beta / norm(beta,"2")
Xquta = X

h = 0.5 #bandwidth
g = log_Y
difference = 1
eps = 1e-8
k = 0

while (difference>eps & k < 1000){
  
  Floor <- sqrt(.Machine$double.eps)
  W = rep(0,N)
  for(i in 1:N){
    if(Delta[i] == 1) W[i] = log_Y[i]
    if(Delta[i] == 0) {
      a = log_Y[i] - g[i]
      W[i] = g[i] + dnorm(a,0,1)/max((1-pnorm(a, 0, 1)), Floor)
    }
  }
  B = W
  CVh = function(param){
    h = param
    h.comp = function(h){
      index = Xquta %*% beta
      W = as.matrix(data.frame(B,1))
      K.sum = npksum(txdat = index, tydat = W, weights = W,bws = h, ckertype = "epanechnikov",leave.one.out = T)$ksum
      g_new = K.sum[1,2,]/NZD(K.sum[2,2,])
      # g = g_hat(beta,X,h,W)
      Floor = sqrt(.Machine$double.eps)
      g_new[which(g_new<Floor)] <- Floor
      g_new[which(g_new>1-Floor)] <- 1 - Floor
      contrib = rep(NA,N)
      for(i in 1:N){
        if(Delta[i]==1) contrib[i] <- (log_Y[i] - g[i])^2
        if(Delta[i]==0){
          a = log_Y[i] - g[i]
          contrib[i] = 1 + (a + 2*g[i] - 2*g_new[i]) * dnorm(a,0,1)/max((1-pnorm(a, 0, 1)), Floor) + (g[i]-g_new[i])^2
        }
      }
      return(sum(contrib))
    }
    if(h>0){return(h.comp(h))}else{return(sqrt(.Machine$double.xmax))}
  }

  h.CV = optim(par = h, CVh, method = "Brent",lower = 0 , upper = 3, control=list(fnscale=1))
  # h.CV = optimize(CVh,c(0,5))
  h = h.CV$par
  bdwth = h
  
  l_I=function(beta){
    Floor <- sqrt(.Machine$double.eps)
    index <- Xquta %*% beta
    W = as.matrix(data.frame(B,1))
    K.sum = npksum(txdat=index, tydat=W,weights=W,bws=h,ckertype="epanechnikov",leave.one.out = T)$ksum
    g_new = K.sum[1,2,]/NZD(K.sum[2,2,])
    # g = g_hat(beta,X,h,B)
    g_new[which(g_new<Floor)] <- Floor
    g_new[which(g_new>1-Floor)] <- 1 - Floor
    contrib <- rep(NA,N)
    for(i in 1:N){
      if(Delta[i]==1) contrib[i] <- (log_Y[i] - g[i])^2
      if(Delta[i]==0){
        a = log_Y[i] - g[i]
        contrib[i] = 1 + (a + 2*g[i] - 2*g_new[i]) * dnorm(a,0,1)/max((1-pnorm(a, 0, 1)), Floor) + (g[i]-g_new[i])^2
      } 
    }
    return(mean(contrib))
  }
  
  para_est = optim(par = beta, fn = l_I, method = "BFGS",control = list(fnscale=1))
  beta_est = para_est$par / norm(para_est$par, "2")
  
  index = Xquta %*% beta_est
  W = as.matrix(data.frame(B,1))
  K.sum = npksum(txdat=index, tydat=W,weights=W,bws= bdwth,ckertype="epanechnikov",leave.one.out = T)$ksum
  g.hat = K.sum[1,2,]/NZD(K.sum[2,2,])
  # g.hat = g_hat(para_est$par,X,h,W)
  
  difference = sum(beta_est - beta)^2 + sum(g - g.hat)^2
  # difference = sum(g - g.hat)^2
  beta = beta_est
  g = g.hat
  k = k+1
}

data.frame(log_Y,g+error)[1:20,]

data.frame(beta_true,beta_est)
