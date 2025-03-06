library(compiler)
library(Rfast)
library(SIS)

###############################################
# Smaller Functions
###############################################
inv_logit <- function(z){
  return(1/(1+exp(-z)))
}

nullpower_append <- function(dataset, power){
  iter = nrow(dataset)
  dataset = dataset[,-c(1,2)]
  take = ncol(power)
  dataset = dataset[,c(1:(take-1))]
  nnrow = colSums(abs(dataset) > 1.96 * 1)/iter
  nptest = rbind(c(0,nnrow),power)
  return(nptest)
}

z_intercept <- function(x, y, b0 = 0){
  n = ncol(x)
  p = nrow(x)
  beta0 = rep(b0, p)
  z_mat = matrix(data = NA, nrow = p, ncol = (n*(n-1)/2))
  index = Rfast::comb_n(n, 2)
  
  for(i in 1:ncol(index)){
      z_mat[,i] = (x[,index[1,i]] - x[,index[2,i]]) * (y[index[1,i]] - y[index[2,i]] - sum(beta0 *  (x[,index[1,i]] - x[,index[2,i]])))
  }
  
  return(z_mat)
}
z_intercept <- compiler::cmpfun(z_intercept)
###############################################
# Bigger Functions
###############################################

# Simulates the necessary data with the specified error
simulate_data <- function(n, c, rho = 0.5, beta = 0, delta = 0, signal = "full", epsilon = "normal", sigma_type = "distance", b0 = 0, xgen = "normal"){
  
  # generating all things that are not x (beta, epsilon)
  if(signal == "full"){beta_vec_samp = c(rep(delta,beta), rep(0, n*c - beta))}
  if(signal == "runif"){beta_vec_samp = c(delta * runif(beta, min = 0, max = 1), rep(0, n*c - beta))}
  #beta_vec_samp = c(rep(delta,beta), rep(0, n*c - beta))
  beta_vec = sample(beta_vec_samp, replace = FALSE)
  beta_vec = beta_vec - rep(0,n*c) # this is to represent the conversion to beta* in line with the paper
  
  # epsilon itself can either be the normal, gamma, or t
  if(epsilon == "normal"){eps_vec = rnorm(n, mean = 0, sd = 1)}
  if(epsilon == "t"){eps_vec = ((4/6)^(1/2))*rt(n, df = 6)}
  
  # generating x, y and z storage variables
  x_mat = matrix(data = 0, nrow = n*c, ncol = n)
  y_vec = rep(0, n)
  z_mat = matrix(data = NA, nrow = n*c, ncol = n)
  
  ## Two different X's depending on the covariance structure
  # This one is rho^{|i-j|}
  if(sigma_type == "distance"){
    for(i in 1:n){
      if(xgen == "normal"){evx = rnorm(n*c, mean = 0, sd = 1)}
      if(xgen == "t"){evx = ((4/6)^(1/2))*rt(n*c, df = 6)}
      for(j in 1:(n*c)){
        if(j == 1){x_mat[j,i] = evx[j]}
        else{
          x_mat[j,i] =  rho * x_mat[j-1, i] + sqrt(1 - rho^2) * evx[j]
        }
      }
      y_vec[i] = sum(x_mat[,i] * beta_vec) + eps_vec[i] + b0
      # NOTE: the z calculated here assumes we know the true means of x and y
      z_mat[,i] = x_mat[,i] * sum(x_mat[,i] * beta_vec) + x_mat[,i] * eps_vec[i]
    }
  }
  # This one is 1 if i = j and rho if i != j
  if(sigma_type == "unit"){
    for(i in 1:n){
      
      if(xgen == "normal"){
        evx = rnorm(n*c, mean = 0, sd = 1)
        z1 = rnorm(1, mean = 0, sd = 1)
      }
      if(xgen == "t"){
        evx = ((4/6)^(1/2))*rt(n*c, df = 6)
        z1 = ((4/6)^(1/2))*rt(1, df = 6)
      }
      for(j in 1:(n*c)){
        x_mat[j,i] =  sqrt(rho) * z1 + sqrt(1 - rho) * evx[j]
      }
      y_vec[i] = sum(x_mat[,i] * beta_vec) + eps_vec[i] + b0
      # NOTE: the z calculated here assumes we know the true means of x and y
      z_mat[,i] = x_mat[,i] * sum(x_mat[,i] * beta_vec) + x_mat[,i] * eps_vec[i]
    }
  }
  
  return(list(z_mat, x_mat, y_vec))
}
simulate_data <- compiler::cmpfun(simulate_data)
###############################################

# Simulates the data for logistic regression
simulate_data_logistic <- function(n, c, rho = 0.5, beta = 0, delta = 0, signal = "full", b0 = 0, xgen = "normal", sigma_type = "distance"){
  
  # generating all things that are not x (beta, epsilon)
  if(signal == "full"){beta_vec_samp = c(rep(delta,beta), rep(0, n*c - beta))}
  if(signal == "runif"){beta_vec_samp = c(delta * runif(beta, min = 0, max = 1), rep(0, n*c - beta))}
  #beta_vec_samp = c(rep(delta,beta), rep(0, n*c - beta))
  beta_vec = sample(beta_vec_samp, replace = FALSE)
  beta_vec = beta_vec - rep(0,n*c) # this is to represent the conversion to beta* in line with the paper
  
  # generating an x, y, and z
  x_mat = matrix(data = 0, nrow = n*c, ncol = n)
  y_vec = rep(0, n)
  z_mat = matrix(data = NA, nrow = n*c, ncol = n)
  
  ## Two different X's depending on the covariance structure
  # This one is rho^{|i-j|}
  if(sigma_type == "distance"){
    for(i in 1:n){
      if(xgen == "normal"){evx = rnorm(n*c, mean = 0, sd = 1)}
      if(xgen == "t"){evx = ((4/6)^(1/2))*rt(n*c, df = 6)}
      for(j in 1:(n*c)){
        if(j == 1){x_mat[j,i] = evx[j]}
        else{
          x_mat[j,i] =  rho * x_mat[j-1, i] + sqrt(1 - rho^2) * evx[j]
        }
      }
      y_vec[i] = rbinom(n = 1, size = 1,inv_logit(sum(x_mat[,i] * beta_vec)))
      z_mat[,i] = x_mat[,i] * (y_vec[i] - 0.5)
    }
  }
  # This one is 1 if i = j and rho if i != j
  if(sigma_type == "unit"){
    for(i in 1:n){
      if(xgen == "normal"){
        evx = rnorm(n*c, mean = 0, sd = 1)
        z1 = rnorm(1, mean = 0, sd = 1)
      }
      if(xgen == "t"){
        evx = ((4/6)^(1/2))*rt(n*c, df = 6)
        z1 = ((4/6)^(1/2))*rt(1, df = 6)
      }
      for(j in 1:(n*c)){
        x_mat[j,i] =  sqrt(rho) * z1 + sqrt(1 - rho) * evx[j]
      }
      y_vec[i] = rbinom(n = 1, size = 1,inv_logit(sum(x_mat[,i] * beta_vec)))
      z_mat[,i] = x_mat[,i] * (y_vec[i] - 0.5)
    }
  }
  
  return(list(z_mat, x_mat, y_vec))
}
simulate_data_logistic <- compiler::cmpfun(simulate_data_logistic)
###############################################


# Simulates the data for poisson regression
simulate_data_poisson <- function(n, c, rho = 0.5, beta = 0, delta = 0, signal = "full", b0 = 0, xgen = "normal", sigma_type = "distance"){
  
  # generating all things that are not x (beta, epsilon)
  if(signal == "full"){beta_vec_samp = c(rep(delta,beta), rep(0, n*c - beta))}
  if(signal == "runif"){beta_vec_samp = c(delta * runif(beta, min = 0, max = 1), rep(0, n*c - beta))}
  #beta_vec_samp = c(rep(delta,beta), rep(0, n*c - beta))
  beta_vec = sample(beta_vec_samp, replace = FALSE)
  beta_vec = beta_vec - rep(0,n*c) # this is to represent the conversion to beta* in line with the paper
  
  # generating an x and z
  x_mat = matrix(data = 0, nrow = n*c, ncol = n)
  y_vec = rep(0, n)
  z_mat = matrix(data = NA, nrow = n*c, ncol = n)
  
  ## Two different X's depending on the covariance structure
  # This one is rho^{|i-j|}
  if(sigma_type == "distance"){
    for(i in 1:n){
      if(xgen == "normal"){evx = rnorm(n*c, mean = 0, sd = 1)}
      if(xgen == "t"){evx = ((4/6)^(1/2))*rt(n*c, df = 6)}
      for(j in 1:(n*c)){
        if(j == 1){x_mat[j,i] = evx[j]}
        else{
          x_mat[j,i] =  rho * x_mat[j-1, i] + sqrt(1 - rho^2) * evx[j]
        }
      }
      y_vec[i] = rpois(n = 1, lambda = exp(sum(x_mat[,i] * beta_vec)))
      z_mat[,i] = x_mat[,i] * (y_vec[i] - 1)
    }
  }
  # This one is 1 if i = j and rho if i != j
  if(sigma_type == "unit"){
    for(i in 1:n){
      if(xgen == "normal"){
        evx = rnorm(n*c, mean = 0, sd = 1)
        z1 = rnorm(1, mean = 0, sd = 1)
      }
      if(xgen == "t"){
        evx = ((4/6)^(1/2))*rt(n*c, df = 6)
        z1 = ((4/6)^(1/2))*rt(1, df = 6)
      }
      for(j in 1:(n*c)){
        x_mat[j,i] =  sqrt(rho) * z1 + sqrt(1 - rho) * evx[j]
      }
      y_vec[i] = rpois(n = 1, lambda = exp(sum(x_mat[,i] * beta_vec)))
      z_mat[,i] = x_mat[,i] * (y_vec[i] - 1)
    }
  }
  
  return(list(z_mat, x_mat, y_vec))
}
simulate_data_poisson <- compiler::cmpfun(simulate_data_poisson)
###############################################

# Calculates both T and T0, depending on the script
calculate_statistic <- function(z, l, k, alpha, Sigma = NULL, T0 = FALSE, symmetric = FALSE){
  # elements needed for a.n and xi.n
  if(symmetric == FALSE){
    n = ncol(z)
    z.bar = rowSums(z)/n
  }
  if(symmetric == TRUE){
    n = (1 + sqrt(1 + 4 * 2 *ncol(z)))/2
    z.bar = Rfast::rowsums(z) * (2/(n*(n-1)))
  }
  
  
  # a.n and xi.n, needed for the W
  a.n = l/(norm(z.bar, type = "2")^2 + k * abs(t(alpha) %*% (z.bar))^2)^(1/2)
  xi.n = (n + 2)/(1 + a.n)
  
  # W, the Sigma.hat, and the trace.hats, needed for the T
  W = -2*(n*log(1 + 1/n*(1-(1+(n/(n+2))*xi.n^2)^(1/2))) + 
            log((1/2) + xi.n/2 + (1/2)*(1 + (n/(n+2)) * xi.n^2)^(1/2)) + 
            log((1/2) - xi.n/2 + (1/2)*(1 + (n/(n + 2)) * xi.n^2)^(1/2)))

  Sigma.hat = Rfast::cova(t(z))
  
  ### Calculating the Trace of Sigma.hat and Sigma.hat^2
  # Gets the matrices that represent x_i^T x_j where i != j
  mlist =  apply(z, FUN = function(x) colSums(x * z), MARGIN = 2)
  diag(mlist) = 0
  
  ## Trace of Sigma.hat
  tr.Sigma.hat = sum(z^2)/n - (1/(n*(n-1)))*sum(mlist)
  
  ## Trace of Sigma.hat^2
  # The ijk case (x_i^T x_j x_j^T x_k, where i,j, and k are all different) and ijkl
  ijk = rep(NA, n) # Empty vector used for storage
  ijkl.store = rep(NA,n)
  sm = sum(mlist)
  cs = Rfast::colsums(mlist)
  for(i in 1:n){
    ijk[i] = sum(mlist[i,] %*% mlist[,-i])
    ijkl.store[i] = sum(mlist[i,] * sm) - sum(2 * mlist[i,] * cs) - sum(2*mlist[i,] * cs[i]) + sum(2*mlist[i,]^2)
  }
  # Summations to get the final values for the HEL equation
  sijk = sum(ijk)
  sijkl = sum(ijkl.store)
  ij = sum(mlist^2)
  
  tr.Sigma.hat2 = (1/(n * (n-1))) * ij - 
    (2/(n*(n-1)*(n-2))) * sijk + 
    (1/(n*(n-1)*(n-2)*(n-3))) * sijkl
  
  # Getting the Omega.hat for the statistic
  tr.omega.hat = tr.Sigma.hat +  k * t(alpha) %*% Sigma.hat %*% alpha
  tr.omega2.hat =  tr.Sigma.hat2 + 
    2*k * t(alpha) %*% Sigma.hat %*% Sigma.hat %*% alpha + 
    k^2 * (t(alpha) %*% Sigma.hat %*% alpha)^2
  
  # T statistic
  T = (2 * tr.omega2.hat)^(-1/2) * ((2 * n * l^2)/((n+2)^2) * W - tr.omega.hat)
  
  if(T0 == TRUE){
    # Calculating the Trace of the Sigma^2 here
    holder = rep(0, nrow(Sigma))
    for(i in 1:nrow(Sigma)){
      holder[i] = sum(Sigma[i,] * Sigma[,i])
    }
    tr.Sigma2 = sum(holder)
    
    tr.omega = sum(diag(Sigma)) +  k * t(alpha) %*% Sigma %*% alpha
    tr.omega2 = tr.Sigma2 + 
      2*k + t(alpha) %*% Sigma %*% Sigma %*% alpha + 
      k^2 * (t(alpha) %*% Sigma %*% alpha)^2
    
    # T.0 statistic 
    T.0 = (2 * tr.omega2)^(-1/2) * ((2 * n * l^2)/((n+2)^2) * W - tr.omega)
    return(c(T, T.0))
  }
  
  return(T)
}
calculate_statistic <- compiler::cmpfun(calculate_statistic)
###############################################

# Calculates both T and T0, depending on the script, based on new method from Changcheng
calculate_statistic_new <- function(z, x, y, kn, alpha, Sigma = NULL, T0 = FALSE, centerx = TRUE){
  # elements needed for W
  n = ncol(z)
  z.bar = rowSums(z)/n
  y.bar = mean(y)
  x.bar = rowSums(x)/n
  
  # this section is for using z.unknown = (y-y.bar)*(x - x.bar)
  y.adjust = y - y.bar
  if(centerx == TRUE){x.adjust = x - x.bar}
  if(centerx == FALSE){x.adjust = x}
  z.unknown = matrix(data = 0, nrow = nrow(x.adjust), ncol = ncol(x.adjust))
  for(iterate in 1:ncol(x.adjust)){
    z.unknown[,iterate] = y.adjust[iterate] * x.adjust[,iterate]
  }
  
  # calculation of the z * z^T or equivalent
  wstore = 0
  wstore.unknown = 0
  for(iterate in 2:n){
    # known mean of x and y
    wstore = wstore + sum(z[,iterate] * z[,1:(iterate-1)]) + kn*sum(alpha*z[,1:(iterate-1)])*sum(z[,iterate]*alpha)
    
    # unknown means of x and y
    wstore.unknown = wstore.unknown + sum(z.unknown[,iterate] * z.unknown[,1:(iterate-1)]) + kn*sum(alpha*z.unknown[,1:(iterate-1)])*sum(z.unknown[,iterate]*alpha)
  }
  
  # This section is for the full method that has fewer assumptions but is more difficult to calculate, however, I am commenting this out to use the simpler method to get simulation results
  #   zsumstore = rep(0,n-1)
  #   wstore.unknown = rep(0,n-1)
  #   y_ij_sumless = rep(0,n)
  #   
  #   for(j in 2:n){
  #     y_ij_sumless[j] =  (2/(n*(n-1)))*sum(y[1:(j-1)]*y[j])
  #     deltay.vec = y[1:(j-1)] * y[j] - (y[1:(j-1)] * (n * y.bar - y[1:(j-1)]))/(n-1) - (y[j]*(n*y.bar - y[j]))/(n-1) + y_ij_sumless[j]
  #     # this basically is to replace the trace sums
  #     trdeltax.vec = t(x[,j]) %*% x[,1:(j-1)] - (1/(n-1)) * colSums(as.matrix((n*x.bar - x[,1:(j-1)])*x[,1:(j-1)])) - as.numeric((1/(n-1)) * (t(n*x.bar - x[,j]) %*% x[,j])) +  sum(diag((2/(n*(n-1)))*Rfast::rowsums(as.matrix(x[,1:(j-1)])) %*% t(x[,j])))
  #     
  #     # this is the second part of the W_n
  #     wxpart2.vec = kn * (colSums(as.matrix(alpha * x[,1:(j-1)])) * sum(alpha * x[,j]) - (1/(n-1))*colSums(as.matrix(alpha * x[,1:(j-1)]))*colSums(as.matrix(alpha * (n*x.bar -x[,1:(j-1)]))) - (1/(n-1)) * sum(alpha * x[,j]) * sum(alpha * (n*x.bar - x[,j])) + t(alpha) %*% ((2/(n*(n-1)))*Rfast::rowsums(as.matrix(x[,1:(j-1)])) %*% t(x[,j])) %*% alpha)
  #     
  #     # adds them all up for each
  #     wstore.unknown[j-1] = (deltay.vec * (trdeltax.vec + wxpart2.vec))[j-1]
  #     
  #     # legacy code for validation purposes, do not use
  # #    for(i in 1:(j-1)){
  # #      deltax = x[,i] %*% t(x[,j]) - (x[,i] %*% t(n*x.bar - x[,i]))/(n-1) - (x[,j] %*% t(n*x.bar-x[,j]))/(n-1) + (2/(n*(n-1)))*Rfast::rowsums(as.matrix(x[,1:(j-1)])) %*% t(x[,j])
  # #      deltaz = deltay.vec[i] * deltax
  # #      zsumstore[j-1] = sum(diag(deltaz)) + kn * t(alpha) %*% deltaz %*% alpha
  # #    }
  #   }
  
  # The W's
  W.pop = 2/(n*(n-1))*wstore
  W.samp = 2/(n*(n-1))*sum(wstore.unknown)
  
  ### Parts for calculating the variance under the known mean case
  ## \hat{\Sigma_z}, sample covariance of z
  sigmaz.hat = cov(t(z))
  
  ## \hat{\Sigma_z}^2
  diagmats = z %*% t(z) # gets the sum of the z_i z_i^T matrices 
  zizjt_sum_ineqj = Rfast::rowsums(z) %*% t(Rfast::rowsums(z)) - diagmats # gets the sum of z_i z_j^T for all i != j
  
  # gets \sum_{i \neq j}(z_j^T z_i)^2
  zjtzi_sum =  apply(z, FUN = function(a) colSums(a * z), MARGIN = 2) # gets the matrices that represent z_i^T z_j for all i and j, including i = j
  diag(zjtzi_sum) = 0 # sets all i = j to 0
  zjtzi_sum_ineqj_sq = sum(zjtzi_sum^2) # final result
  
  # gets \sum_{i \neq j} \alpha^T (z_i z_j^T)^2 \alpha
  all_alphatz = colSums(alpha * z)
  alphatzizjt2alpha_sum_ineqj = sum(all_alphatz * t(all_alphatz * zjtzi_sum))
  
  ## Final \hat{\tr(\Omega^2)} to be used in variance
  tr.omega.hat.sq.known = 1/(n*(n-1)) * zjtzi_sum_ineqj_sq + sum(diag((1/(n*(n-1)))^2 * zizjt_sum_ineqj %*% zizjt_sum_ineqj)) + (2*kn)/(n*(n-1)) * alphatzizjt2alpha_sum_ineqj -  (2*kn)/(n*(n-1)) * t(alpha) %*% ((1/(n*(n-1)))^2 * zizjt_sum_ineqj %*% zizjt_sum_ineqj) %*% alpha + kn^2 * (t(alpha) %*% sigmaz.hat %*% alpha)^2
  
  ### Parts for calculating the variance under the unknown mean case
  ## \hat{\Sigma_z}^2
  diagmats.unknown = z.unknown %*% t(z.unknown) # gets the sum of the z_i z_i^T matrices 
  zizjt_sum_ineqj.unknown = Rfast::rowsums(z.unknown) %*% t(Rfast::rowsums(z.unknown)) - diagmats.unknown # gets the sum of z_i z_j^T for all i != j
  
  # gets \sum_{i \neq j}(z_j^T z_i)^2
  zjtzi_sum.unknown =  apply(z.unknown, FUN = function(a) colSums(a * z.unknown), MARGIN = 2) # gets the matrices that represent z_i^T z_j for all i and j, including i = j
  diag(zjtzi_sum.unknown) = 0 # sets all i = j to 0
  zjtzi_sum_ineqj_sq.unknown = sum(zjtzi_sum.unknown^2) # final result
  
  # gets \sum_{i \neq j} \alpha^T (z_i z_j^T)^2 \alpha
  all_alphatz.unknown = colSums(alpha * z.unknown)
  alphatzizjt2alpha_sum_ineqj.unknown = sum(all_alphatz.unknown * t(all_alphatz.unknown * zjtzi_sum.unknown))
  
  ## Final \hat{\tr(\Omega^2)} to be used in variance
  tr.omega.hat.sq.unknown = 1/(n*(n-1)) * zjtzi_sum_ineqj_sq.unknown + sum(diag((1/(n*(n-1)))^2 * zizjt_sum_ineqj.unknown %*% zizjt_sum_ineqj.unknown)) + (2*kn)/(n*(n-1)) * alphatzizjt2alpha_sum_ineqj.unknown -  (2*kn)/(n*(n-1)) * t(alpha) %*% ((1/(n*(n-1)))^2 * zizjt_sum_ineqj.unknown %*% zizjt_sum_ineqj.unknown) %*% alpha + kn^2 * (t(alpha) %*% sigmaz.hat %*% alpha)^2
  
  ## This section below I can't get to work just yet so we're using a simpler method for now
  # ## for tr(\hat{\Sigma}_z^2)
  # # this part for the \delta_{ij}(y)
  # y_ij_mats = y %*% t(y) # gets all the y_i y_j
  # y_subtract = y*(n*(y.bar) - y)/(n-1) # creates the y_i(n * y.bar - y_i)/(n-1) stuff, which is the same for the j's
  # y_ij_mats = y_ij_mats - y_subtract
  # y_ij_mats = sweep(y_ij_mats, 2, y_subtract)
  # diag(y_ij_mats) = 0 # removes all y_i y_i
  # for(iter in 1:n){
  #   if(iter == 1){y_ij_mats[iter,] = y_ij_mats[iter,] + y_ij_sumless}
  #   else{
  #     y_ij_mats[iter,] = y_ij_mats[iter,] + c(y_ij_sumless[2:iter],y_ij_sumless[-(2:iter)])
  #   }
  # }
  # 
  # # this part for the \delta_{ij}(x)
  # # subtraction part
  # minus_ij = sum(y_ij_mats)*(n*(n-1)/2)*(1/(n-1))*(n*Rfast::rowsums(x) %*% t(x.bar) + x %*% t(x))
  # # uh the sum part????
  # weirdstore = matrix(0, nrow = nrow(x), ncol = nrow(x))
  # for(j in 2:n){
  #   weirdstore = weirdstore + (2/(n*(n-1)))*Rfast::rowsums(as.matrix(y_ij_mats[1:(j-1),j] * x[,1:(j-1)])) %*% t(x[,j])
  # }
  # 
  # diagmats_x = x %*% t(x) # gets the sum of the x_i x_i^T matrices 
  # xixjt_sum_ineqj = Rfast::rowsums(x %*% y_ij_mats) %*% t(Rfast::rowsums(x)) - diagmats_x # gets the sum of z_i z_j^T for all i != j
  # deltaz.mat = xixjt_sum_ineqj - minus_ij + weirdstore
  # 
  # ## Final \hat{\tr(\Omega^2)} to be used in variance
  # tr.omega.hat.sq.known = 1/(n*(n-1)) * zjtzi_sum_ineqj_sq + sum(diag((1/(n*(n-1)))^2 * deltaz.mat %*% deltaz.mat)) + (2*kn)/(n*(n-1)) * alphatzizjt2alpha_sum_ineqj - t(alpha) %*% ((1/(n*(n-1)))^2 * deltaz.mat %*% deltaz.mat) %*% alpha + kn^2 * (t(alpha) %*% sigmaz.hat %*% alpha)^2
  
  # T statistic for known means
  T.pop = (2 * tr.omega.hat.sq.known)^(-1/2) * n * W.pop
  T.samp = (2 * tr.omega.hat.sq.unknown)^(-1/2) * n * W.samp
  
  # If we pass the true covariance matrix for z
  if(T0 == TRUE){
    # Getting Sigma_z^2
    sigmaz.sq = Sigma %*% Sigma
    tr.omega.sq.true = sum(diag(sigmaz.sq)) + 2*kn*t(alpha) %*% sigmaz.sq %*% alpha + kn^2 *(t(alpha) %*% Sigma %*% alpha)^2
    T0.pop = (2 * tr.omega.sq.true)^(-1/2) * n * W.pop
    T0.samp = (2 * tr.omega.sq.true)^(-1/2) * n * W.samp
    return(c(T.pop, T0.pop, T.samp, T0.samp))
  }
  
  return(c(T.pop, T.samp))
}
calculate_statistic_new <- compiler::cmpfun(calculate_statistic_new)
###############################################

# Zhong and Chen 2011 comparison on both T and T0, depending on the script
# NOTE: this is based on the symmetrized version of T_np found under section 4 of main results (p263 or p5 of the pdf)
zhongchen2011 <- function(X, y, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE){
  # elements needed for a statistic
  n = ncol(X)
  p = nrow(X)
  small.sigma = var(y)
  beta0 = c(rep(delta,beta), rep(0, p - beta)) 
  
  ### Calculating the Trace of Sigma.hat^2
  # Gets the matrices that represent x_i^T x_j where i != j
  mlist =  apply(X, FUN = function(x) colSums(x * X), MARGIN = 2)
  diag(mlist) = 0
  
  ## Trace of Sigma.hat^2 and T_np at the same time
  # ijkl cases and storage stuff for trace
  ijk = rep(0, n) # Empty vector used for storage
  ijkl.store = rep(0,n)
  sm = sum(mlist)
  cs = Rfast::colsums(mlist)
  
  # storage for T_np
  diffy = matrix(data = 0, nrow = n, ncol = n)
  dxy = matrix(data = 0, nrow = p, ncol = n)
  dxy2 = matrix(data = 0, nrow = p, ncol = n)
  horizontal_all_arrays_sum = rep(0, p)
  vertical_all_arrays_sum = matrix(data = 0, nrow = p, ncol = n)
  horizontal_within_arrays = matrix(data = 0, nrow = p, ncol = n)
  
  for(i in 1:n){
    ## Trace of Sigma stuff
    ijk[i] = sum(mlist[i,] %*% mlist[,-i])
    ijkl.store[i] = sum(mlist[i,] * sm) - sum(2 * mlist[i,] * cs) - sum(2*mlist[i,] * cs[i]) + sum(2*mlist[i,]^2)
    
    ## T_np stuff
    ## Getting every (X_ia - X_ib) and (y_ia - y_ib - (X_ia - X_ib)^T * Beta0) combination
    dxy2 = X[,i] - X
    diffy[i,] = y[i] - y - Rfast::colsums(dxy2 * beta0)
    
    # one weird fix
    dxy2[ ,1:i] = 0 
    diffy[i,1:i] = 0 
    
    ## Combines the two together
    dxy2 = t(diffy[i,] * t(dxy2))
    
    ## all the storage stuff
    horizontal_within_arrays[,i] = Rfast::rowsums(dxy2)
    vertical_all_arrays_sum = vertical_all_arrays_sum + dxy2
    
    # this replaces the sum(dxy^2) below from when it was an array
    dxy = dxy + dxy2^2
  }
  
  
  # Summations to get the final values for the variance equation
  sijk = sum(ijk)
  sijkl = sum(ijkl.store)
  ij = sum(mlist^2)
  
  tr.Sigma.hat2 = (1/(n * (n-1))) * ij - 
    (2/(n*(n-1)*(n-2))) * sijk + 
    (1/(n*(n-1)*(n-2)*(n-3))) * sijkl
  
  # gets the product of all the (X_ia - X_ib)^T (X_ic - X_id)(y_ia - y_ib - (X_ia - X_ib)^T * Beta0)(y_ic - y_id - (X_ic - X_id)^T * Beta0)
  horizontal_all_arrays_sum = Rfast::rowsums(horizontal_within_arrays)
  
  newresult = sum(horizontal_all_arrays_sum^2) - sum(vertical_all_arrays_sum^2) -  sum(horizontal_within_arrays^2) - 2*sum(vertical_all_arrays_sum *horizontal_within_arrays) + sum(dxy)
  
  # Zhong and Chen statistic
  T = (1/4) * (1/3) * (4*3*2)/(n*(n-1)*(n-2)*(n-3)) * newresult/2 #divide by two because the newresult construction has (x_1 - x_2)(x_3- x_4) and (x_3 - x_4)(x_1 - x_2) both, so you need to halve it
  # varT = (2/(n*(n-1)))*tr.Sigma.hat2 * small.sigma^2
  normalized.T = (n*T)/(sqrt(2*tr.Sigma.hat2) * small.sigma)
  
  if(T0 == TRUE){
    # Calculating the Trace of the Sigma^2 here
    holder = rep(0, nrow(Sigma))
    for(i in 1:nrow(Sigma)){
      holder[i] = sum(Sigma[i,] * Sigma[,i])
    }
    tr.Sigma2 = sum(holder)
    
    # Zhong and Chen statistic with true sigma 
    normalized.T0 = (n*T)/(sqrt(2*tr.Sigma2) * small.sig)
    return(c(T, normalized.T, normalized.T0))
  }
  
  return(c(T,normalized.T))
}
zhongchen2011 <- compiler::cmpfun(zhongchen2011)
###############################################

# Cui et al 2018 comparison on both T and T0, depending on the script
cui2018 <- function(X, y, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE){
  # elements needed for a statistic
  n = ncol(X)
  small.sigma = var(y)
  beta0 = c(rep(delta,beta), rep(0, n*c - beta)) 
  
  ### RCV
  # randomly split the data in half
  split = sample(1:n, round(n/2), replace = FALSE)
  
  X1 = t(X[,split])
  X2 = t(X[,-split])
  y1 = y[split]
  y2 = y[-split]
  
  #SIS
  sis1_m2 = SIS::SIS(X2, y2, family = "gaussian", penalty = "lasso", iter = FALSE)
  sis1x = X1[,sis1_m2$sis.ix0]
  s1 = (t(y1) %*% (diag((n - round(n/2))) - (sis1x) %*% solve(t(sis1x) %*% (sis1x)) %*% t(sis1x)) %*% (y1))/(round(n/2) - length(sis1_m2$sis.ix0))

  sis2_m1 = SIS::SIS(X1, y1, family = "gaussian", penalty  = "lasso", iter = FALSE)
  sis2x = X1[,sis2_m1$sis.ix0]
  s2 = (t(y2) %*% (diag(round(n/2)) - (sis2x) %*% solve(t(sis2x) %*% (sis2x)) %*% t(sis2x)) %*% (y2))/(n/2 - length(sis2_m1$sis.ix0))
  
  sig.rcv = (s1 + s2)/2
  
  ### Calculating the Trace of Sigma.hat^2
  # Gets the matrices that represent x_i^T x_j where i != j
  mlist =  apply(X, FUN = function(x) colSums(x * X), MARGIN = 2)
  diag(mlist) = 0
  
  ## Trace of Sigma.hat^2 and T_np at the same time
  # ijkl cases and storage stuff for trace
  ijk = rep(0, n) # Empty vector used for storage
  ijkl.store = rep(0,n)
  sm = sum(mlist)
  cs = Rfast::colsums(mlist)
  
  # T_np storage stuff and precalculations
  xbar = Rfast::rowsums(X)/n
  ybar = mean(y)
  doublesum = 0
  
  for(i in 1:n){
    ## Trace of Sigma stuff
    ijk[i] = sum(mlist[i,] %*% mlist[,-i])
    ijkl.store[i] = sum(mlist[i,] * sm) - sum(2 * mlist[i,] * cs) - sum(2*mlist[i,] * cs[i]) + sum(2*mlist[i,]^2)
    
    ## T_np stuff
    if(i == 1){
      doublesum = doublesum + 0
    } else{
      for(j in 1:(i-1)){
        deltax = sum((X[,i] - xbar) * (X[,j] - xbar)) + sum((X[,i] - X[,j])^2)/(2*n)
        deltay = (y[i] - ybar) * (y[j] - ybar) + ((y[i] - y[j])^2)/(2*n)
        doublesum = doublesum + deltax * deltay
      } 
    }
  }
  
  # Summations to get the final values for the variance equation
  sijk = sum(ijk)
  sijkl = sum(ijkl.store)
  ij = sum(mlist^2)
  
  tr.Sigma.hat2 = (1/(n * (n-1))) * ij - 
    (2/(n*(n-1)*(n-2))) * sijk + 
    (1/(n*(n-1)*(n-2)*(n-3))) * sijkl
  
  # Cui statistic
  T = (1 - 2/n)^(-2) * (2/(n*(n-1))) * doublesum
  normalized.T = (n*T)/(sqrt(2*tr.Sigma.hat2) * small.sigma)
  rcv.T = (n*T)/(sqrt(2*tr.Sigma.hat2) * sig.rcv)
  
  if(T0 == TRUE){
    # Calculating the Trace of the Sigma^2 here
    holder = rep(0, nrow(Sigma))
    for(i in 1:nrow(Sigma)){
      holder[i] = sum(Sigma[i,] * Sigma[,i])
    }
    tr.Sigma2 = sum(holder)
    
    # Zhong and Chen statistic with true sigma 
    normalized.T0 = (n*T)/(sqrt(2*tr.Sigma2) * small.sig)
    rcv.T0 = (n*T)/(sqrt(2*tr.Sigma2) * sig.rcv)
    return(c(T, normalized.T, normalized.T0, rcv.T, rcv.T0))
  }
  return(c(T, normalized.T, rcv.T))
}
cui2018 <- compiler::cmpfun(cui2018)