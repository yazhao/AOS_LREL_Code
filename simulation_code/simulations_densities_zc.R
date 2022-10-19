###############################################
# This script simulates the data according to 
# various null hypotheses (different errors),
# calculates the corresponding statistics T and
# T_0, and then plots their densities versus a
# standard normal distribution
#
# Because the null hypothesis is that
# Beta = Beta_0, the end result is that
# y = epsilon with mean 0 and variance sigma^2
#
# For convenience we've set sigma^2 = 1
###############################################

###############################################
# Libraries, function sourcing, parameters
###############################################
source("simulation_code/simulations_functions.R")
dir = "results/distance/"

n = c(200, 400, 800)
c = 1.5
iter = 1000
cov_rho = 0.5

###############################################

###############################################
# Errors are Normal N(0,1)
###############################################
# Storage
normal.null = matrix(data = NA, nrow = iter, ncol = 2*length(n)+1)
normal.null[,1] = 1:iter

for(i in 1:length(n)){
  # The true Sigma, only used for these densities
  Sigma = matrix(data = 0, nrow = n[i] * c, ncol = n[i] * c)
  Sigma = cov_rho^abs(row(Sigma) - col(Sigma))
  Sigma = 1 * Sigma
  
  for(k in 1:iter){
    
    # Setting the seed 
    set.seed(k)

    # generating the simulated data
    results = simulate_data(n = n[i], c = c, rho = cov_rho, epsilon = "normal")
    x.mat = results[[2]]
    y.vec = results[[3]]

    # calculating the T statistics
    Tstats = zhongchen2011(X = x.mat, y = y.vec, beta = 0, delta = 0, Sigma = Sigma, small.sig = 1, T0 = TRUE)
    # store in matrix
    normal.null[k,i+1] = Tstats[2]
    normal.null[k,i+1+length(n)] = Tstats[3]
  }
}

write.csv(normal.null, file = paste0(dir,"normal_null_zc_rho", cov_rho ,".csv"))

###############################################

###############################################
# Errors are t-distribution with 6 df
###############################################
# Storage
t.null = matrix(data = NA, nrow = iter, ncol = 2*length(n)+1)
t.null[,1] = 1:iter

for(i in 1:length(n)){
  # The true Sigma, only used for these densities
  Sigma = matrix(data = 0, nrow = n[i] * c, ncol = n[i] * c)
  Sigma = cov_rho^abs(row(Sigma) - col(Sigma))
  Sigma = 1 * Sigma
  
  for(k in 1:iter){
    
    # Setting the seed 
    set.seed(k)
    
    # generating the simulated data
    results = simulate_data(n = n[i], c = c, rho = cov_rho, epsilon = "t")
    x.mat = results[[2]]
    y.vec = results[[3]]
    
    # calculating the T statistics
    Tstats = zhongchen2011(X = x.mat, y = y.vec, beta = 0, delta = 0, Sigma = Sigma, small.sig = 1, T0 = TRUE)
    
    # store in matrix
    t.null[k,i+1] = Tstats[2]
    t.null[k,i+1+length(n)] = Tstats[3]
  }
}

write.csv(t.null, file = paste0(dir,"t_null_zc_rho", cov_rho ,".csv"))

###############################################