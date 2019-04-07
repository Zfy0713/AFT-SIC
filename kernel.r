K <- function(u){
  K.res <- rep(NA,length(u))
  for(i in 1:length(u)){
    if(u[i]^2 <= 5)  K.res[i] = 3 * (1-u[i]^2/5)/(4 * sqrt(5))
    else K.res[i] = 0
  }
  return(K.res)
}

g_hat <- function(beta, X, bwh = 0.5, W){
  n = nrow(X)
  g.hat = rep(NA,n)
  for(i in 1:n){
    K.sum = K((X[i,] - X[-i,]) %*% beta / h)
    g.hat[i] = sum(K.sum * W[-i]) / sum(K.sum)
  }
  return(g.hat)
}
