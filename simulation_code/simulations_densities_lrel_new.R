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
cov_rho = 0.25

###############################################

###############################################
# Errors are Normal N(0,1), MEANS KNOWN
###############################################
# Storage
normal.null = matrix(data = NA, nrow = iter, ncol = 2*length(n)+1)
normal.null[,1] = 1:iter
normal.null2 = normal.null

for(i in 1:length(n)){
  # The true Sigma, only used for these densities
  Sigma = matrix(data = 0, nrow = n[i] * c, ncol = n[i] * c)
  Sigma = cov_rho^abs(row(Sigma) - col(Sigma))
  Sigma = 1 * Sigma
  
  # Calculating the k and alpha needed for the T and T_0
  k.n = ((n[i]*c)/log((n[i]*c)))^(1/2)
  alpha = rep(1, (n[i]*c))/sqrt((n[i]*c))
  
  for(k in 1:iter){
    # Setting the seed 
    set.seed(k)

    # generating the simulated data
    results = simulate_data(n = n[i], c = c, rho = cov_rho, epsilon = "normal")
    
    # calculating the T statistics, currently only for known means
    Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = Sigma, T0 = TRUE)
    
    # store in matrix, separately for known and unknown
    normal.null[k,i+1] = Tstats[1]
    normal.null[k,i+1+length(n)] = Tstats[2]
    normal.null2[k,i+1] = Tstats[3]
    normal.null2[k,i+1+length(n)] = Tstats[4]
  }
}

write.csv(normal.null, file = paste0(dir,"knownnew_normal_null_rho", cov_rho ,".csv"))
write.csv(normal.null2, file = paste0(dir,"unknownnew_normal_null_rho", cov_rho ,".csv"))

###############################################

###############################################
# Errors are t-distribution with 6 df, MEANS KNOWN
###############################################
# Storage
t.null = matrix(data = NA, nrow = iter, ncol = 2*length(n)+1)
t.null[,1] = 1:iter
t.null2 = t.null

for(i in 1:length(n)){
  # The true Sigma, only used for these densities
  Sigma = matrix(data = 0, nrow = n[i] * c, ncol = n[i] * c)
  Sigma = cov_rho^abs(row(Sigma) - col(Sigma))
  Sigma = 1 * Sigma
  
  # Calculating the l, k, and alpha needed for the T and T_0
  k.n = ((n[i]*c)/log((n[i]*c)))^(1/2)
  alpha = rep(1, (n[i]*c))/sqrt((n[i]*c))
  
  for(k in 1:iter){
    # Setting the seed 
    set.seed(k)
    
    # generating the simulated data
    results = simulate_data(n = n[i], c = c, rho = cov_rho, epsilon = "t")
    
    # calculating the T statistics
    Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = Sigma, T0 = TRUE)
    
    # store in matrix separately for known and unknown
    t.null[k,i+1] = Tstats[1]
    t.null[k,i+1+length(n)] = Tstats[2]
    t.null2[k,i+1] = Tstats[3]
    t.null2[k,i+1+length(n)] = Tstats[4]
  }
}

write.csv(t.null, file = paste0(dir,"knownnew_t_null_rho", cov_rho ,".csv"))
write.csv(t.null2, file = paste0(dir,"unknownnew_t_null_rho", cov_rho ,".csv"))
##############################################