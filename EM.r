source('GenerateData.R')

# K = function(u){
#   if(u^2 <= 5){
#     k = 3 * (1 - u^2 / 5)/(4 * sqrt(5))
#   }
#   else{ k = 0 }
#   return(k)
# }

beta = matrix(runif(length(mu),-3,3),ncol = 1)
Xquta = X

#beta0 = 0 no intercept
# beta = matrix(runif(length(mu),-3,3),ncol = 1)
# Xquta = X

emmax = 1e6
difference = 1
eps = 1e-8
k = 0
while (difference>eps){
  #E-step
  l_I = function(beta){
    Q = rep(0,N)
    Floor <- sqrt(.Machine$double.eps) 
    for(i in 1:N){
      if(Delta[i] == 1) Q[i] = (log_Y[i] - sum(beta * Xquta[i,]))^2
      if(Delta[i] == 0) {
        # a = log_Y[i] - sum(beta * Xquta[i,])^2
        a = log_Y[i] - sum(beta * Xquta[i,])
        Q[i] = 1 + a * dnorm(a, 0, 1)/max((1-pnorm(a, 0, 1)), Floor)
      }
    }
    
    return(sum(Q))
  }
  
  beta_est = optim(par = beta, l_I, control = list(fnscale = 1))
  
  difference = norm(beta - beta_est$par, type = c("2"))
  beta = beta_est$par
  k = k + 1
}

data.frame(beta,beta_true)

