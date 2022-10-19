###############################################
# This script tests the power settings from
# an earlier script but with different alphas
# and other values
#
# For convenience we've set sigma^2 = 1
###############################################

###############################################
# Libraries, function sourcing, parameters
###############################################
source("simulation_code/simulations_functions.R")
dir = "results/distance/"

n = 400
c = 1.5
delta = c(0, 0.0125, 0.025, 0.0375, 0.05)
delta2 = c(0, 0.1, 0.15, 0.2, 0.3)	
iter = 1000

betashare = 0.25
cov_rho = 0.50

###############################################

###############################################
# Errors are Normal N(0,1)
###############################################

########### Different Alphas ##################
#### Alpha is e_i, where i is the ith element
# of |zbar - mu| that is the largest
# storage for all different deltas and alphas
normal.alpha = matrix(data = NA, nrow = length(delta), ncol = 4)
normal.alpha[,1] = delta
# storing for each of the delta runs
normal.run = matrix(data = 0, nrow = iter, ncol = length(delta))
normal.run2 = normal.run
normal.run3 = normal.run

for(j in 1:length(delta)){
    for(k in 1:iter){
    # Setting the seed and defining the mean and covariance matrix
    set.seed(k)

    # generating the simulated data
    results = simulate_data(n = n, c = c, rho = cov_rho, beta = betashare * n, delta = delta[j], signal = "full", epsilon = "normal")
    z.mat = results[[1]]
    
    # Calculating the k, and alpha needed for the T and T_0
    k.n = ((n*c)/log(n*c))^(1/2)
    
    # regular alpha
    alpha = rep(1, (n * c))/sqrt(n * c)

    # alpha is all zeros except the ith element, which is the maximum of |zbar - mu|, or functionally the max of |zbar|
    z.mat.bar = abs(rowSums(z.mat)/ncol(z.mat))
    alpha.e = rep(0, (n*c))
    alpha.e[which.max(z.mat.bar)] = 1

    # alpha is a uniformly distributed random direction with a norm equal to 1
    alpha.1 = runif((n*c), min = -1, max = 1)
    alpha2 = alpha.1/norm(alpha.1, type = "2")
    
    # calculating the T statistics
    Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha.e, Sigma = NULL, T0 = FALSE)
    Tstats2 = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha2, Sigma = NULL, T0 = FALSE)
    Tstats3 = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
    
    # store in matrix
    normal.run[k,j] = Tstats[1]
    normal.run2[k,j] = Tstats2[1]
    normal.run3[k,j] = Tstats3[1]
    }
  
  # calculating power, which is the proportion of test statistics >1.96 (two-sided)
  normal.alpha[j, 2] = sum((abs(normal.run[,j]) > 1.96)*1)/iter
  normal.alpha[j, 3] = sum((abs(normal.run2[,j]) > 1.96)*1)/iter
  normal.alpha[j, 4] = sum((abs(normal.run3[,j]) > 1.96)*1)/iter
}

write.csv(normal.alpha, file = paste0(dir,"sensitivity_alpha.csv"))

########### Different betashare ##################
# storage for all different delta and betashare values
normal.bs = matrix(data = NA, nrow = length(delta), ncol =  4)
normal.bs[,1] = delta

# storing for each of the delta runs
normal.run = matrix(data = 0, nrow = iter, ncol = length(delta))
normal.run2 = normal.run
normal.run3 = normal.run

for(j in 1:length(delta)){
    for(k in 1:iter){
    # Setting the seed
    set.seed(k)
      
    # betashares
    betashare1 = 0.10
    betashare2 = 0.33

    # generating the simulated data each according to the different betashare
    results1 = simulate_data(n = n, c = c, rho = cov_rho, beta = betashare1 * n, delta = delta[j], signal = "full", epsilon = "normal")
    results2 = simulate_data(n = n, c = c, rho = cov_rho, beta = betashare2 * n, delta = delta[j], signal = "full", epsilon = "normal")
    results3 = simulate_data(n = n, c = c, rho = cov_rho, beta = betashare * n, delta = delta[j], signal = "full", epsilon = "normal")
    
    # Calculating the k and alpha needed for the T and T_0
    k.n = ((n*c)/log(n*c))^(1/2)
    alpha = rep(1, (n*c))/sqrt((n*c))
    
    
    # calculating the T statistics
    Tstats = calculate_statistic_new(z = results1[[1]], x = results1[[2]], y = results1[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
    Tstats2 = calculate_statistic_new(z = results2[[1]], x = results2[[2]], y = results2[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
    Tstats3 = calculate_statistic_new(z = results3[[1]], x = results3[[2]], y = results3[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
    
    # store in matrix
    normal.run[k,j] = Tstats[1]
    normal.run2[k,j] = Tstats2[1]
    normal.run3[k,j] = Tstats3[1]
    }
  
  # calculating power, which is the proportion of test statistics >1.96 (two-sided)
  normal.bs[j, 2] = sum((abs(normal.run[,j]) > 1.96)*1)/iter
  normal.bs[j, 3] = sum((abs(normal.run2[,j]) > 1.96)*1)/iter
  normal.bs[j, 4] = sum((abs(normal.run3[,j]) > 1.96)*1)/iter
}

write.csv(normal.bs, file = paste0(dir,"sensitivity_betashare.csv"))

########### Different k.n ##################
# storage for all different delta and k.n values
normal.kn = matrix(data = NA, nrow = length(delta), ncol = 4)
normal.kn[,1] = delta

# storing for each of the delta runs
normal.run = matrix(data = 0, nrow = iter, ncol = length(delta))
normal.run2 = normal.run
normal.run3 = normal.run

for(j in 1:length(delta)){
    for(k in 1:iter){
    # Setting the seed
    set.seed(k)

    # generating the simulated data
    results = simulate_data(n = n, c = c, rho = cov_rho, beta = betashare * n, delta = delta[j], signal = "full", epsilon = "normal")
    
    # Calculating the k and alpha needed for the T and T_0
    k.n1 = 0.5*((n*c)/log(n*c))^(1/2)
    k.n2 = 2.0*((n*c)/log(n*c))^(1/2)
    k.n3 = ((n*c)/log(n*c))^(1/2)
    alpha = rep(1, (n*c))/sqrt(n*c)

    # calculating the T statistics
    Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n1, alpha = alpha, Sigma = NULL, T0 = FALSE)
    Tstats2 = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n2, alpha = alpha, Sigma = NULL, T0 = FALSE)
    Tstats3 = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n3, alpha = alpha, Sigma = NULL, T0 = FALSE)
    
    # store in matrix
    normal.run[k,j] = Tstats[1]
    normal.run2[k,j] = Tstats2[1]
    normal.run3[k,j] = Tstats3[1]
    }
  
  # calculating power, which is the proportion of test statistics >1.96 (two-sided)
  normal.kn[j, 2] = sum((abs(normal.run[,j]) > 1.96)*1)/iter
  normal.kn[j, 3] = sum((abs(normal.run2[,j]) > 1.96)*1)/iter
  normal.kn[j, 4] = sum((abs(normal.run3[,j]) > 1.96)*1)/iter
}

write.csv(normal.kn, file = paste0(dir,"sensitivity_kn.csv"))

###############################################