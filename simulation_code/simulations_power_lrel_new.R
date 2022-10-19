###############################################
# This script simulates both null and 
# alternative hypothesis data and checks the 
# power (which in the null case is the type I
# error). We use a beta defined as a delta *
# [1,...,1,0,...,0], starting with 10% 1's
#
# For convenience we've set sigma^2 = 1
#
# This script is for the LREL method
###############################################

###############################################
# Libraries, function sourcing, parameters
###############################################
source("simulation_code/simulations_functions.R")
dir = "results/distance/"

n = c(200, 400, 800)
c = 1.5
delp = seq(from = 0.25, to = 1, by = 0.25)
iter = 1000
betashare = 0.25
cov_rho = 0.25
# various max delta values for different cov_rho and full vs runif
mdf = c(0.075, 0.05, 0.035)
md25r = c(0.2, 0.1, 0.05)
md5075r = c(0.1,0.05, 0.035)

# this is where you actually set it for mdr and change it
mdr = md25r

###############################################

###############################################
# Errors are Normal, Full Signal
###############################################
# storage for all different delta values
normal.power.array = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
normal.power = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
normal.power[,1] = delp
# storing for each of the delta runs
normal.run = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results from normal_null
  normal_null <- read.csv(paste0(dir,"knownnew_normal_null_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    
    delta = delp * mdf[i]
    # Calculating the k, and alpha needed for the T and T_0
    k.n = ((n[i] * c)/log(n[i] * c))^(1/2)
    alpha = rep(1, (n[i] * c))/sqrt(n[i] * c)
    
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "full", epsilon = "normal")
        # calculating the T statistics
        Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
        
        # store in matrix
        normal.run[k,j] = Tstats[1]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      normal.power[j, i+1] = sum((abs(normal.run[,j]) > 1.96)*1)/iter
    }
  }
  normal.power.array[, , b] = nullpower_append(normal_null,normal.power)
  write.csv(normal.power.array[, , b], file = paste0(dir, paste0("knownnew_normal_power_rho", cov_rho, "signalfull_beta_", betashare[b],".csv")))
}

###############################################

###############################################
# Errors are Normal, Runif Signal
###############################################
# storage for all different delta values
normal.power.array2 = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
normal.power2 = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
normal.power2[,1] = delp
# storing for each of the delta runs
normal.run2 = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results from normal_null
  normal_null <- read.csv(paste0(dir,"knownnew_normal_null_rho", cov_rho, ".csv"))

  for(i in 1:length(n)){
    delta = delp * mdr[i]
    
    # Calculating the k and alpha needed for the T and T_0
    k.n = ((n[i] * c)/log(n[i] * c))^(1/2)
    alpha = rep(1, (n[i] * c))/sqrt(n[i] * c)
    
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "runif", epsilon = "normal")
        
        # calculating the T statistics
        Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
        
        # store in matrix
        normal.run2[k,j] = Tstats[1]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      normal.power2[j, i+1] = sum((abs(normal.run2[,j]) > 1.96)*1)/iter
    }
  }
  normal.power.array2[, , b] = nullpower_append(normal_null,normal.power2)
  write.csv(normal.power.array2[, , b], file = paste0(dir, paste0("knownnew_normal_power_rho", cov_rho, "signalrunif_beta_", betashare[b],".csv")))
}

###############################################

###############################################
# Errors are t-distribution with df = 6, Full Signal
###############################################
# storage for all different delta values
t.power.array = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
t.power = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
t.power[,1] = delp
# storing for each of the delta runs
t.run = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results from normal_null
  t_null <- read.csv(paste0(dir,"knownnew_t_null_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    delta = delp * mdf[i]
    
    # Calculating the k and alpha needed for the T and T_0
    k.n = ((n[i] * c)/log(n[i] * c))^(1/2)
    alpha = rep(1, (n[i] * c))/sqrt(n[i] * c)
    
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "full", epsilon = "t")
        
        # calculating the T statistics
        Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
        
        # store in matrix
        t.run[k,j] = Tstats[1]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      t.power[j, i+1] = sum((abs(t.run[,j]) > 1.96)*1)/iter
    }
  }
  t.power.array[, , b] = nullpower_append(t_null,t.power)
  write.csv(t.power.array[, , b], file = paste0(dir, paste0("knownnew_t_power_rho", cov_rho, "signalfull_beta_", betashare[b],".csv")))
}

###############################################

###############################################
# Errors are t-distribution with df = 6, Runif Signal
###############################################
# storage for all different delta values
t.power.array2 = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
t.power2 = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
t.power2[,1] = delp
# storing for each of the delta runs
t.run2 = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results from normal_null
  t_null <- read.csv(paste0(dir,"knownnew_t_null_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    delta = delp * mdr[i]
    
    # Calculating the k and alpha needed for the T and T_0
    k.n = ((n[i] * c)/log(n[i] * c))^(1/2)
    alpha = rep(1, (n[i] * c))/sqrt(n[i] * c)
    
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "runif", epsilon = "t")
        
        # calculating the T statistics
        Tstats = calculate_statistic_new(z = results[[1]], x = results[[2]], y = results[[3]], kn = k.n, alpha = alpha, Sigma = NULL, T0 = FALSE)
        
        # store in matrix
        t.run2[k,j] = Tstats[1]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      t.power2[j, i+1] = sum((abs(t.run2[,j]) > 1.96)*1)/iter
    }
  }
  t.power.array2[, , b] = nullpower_append(t_null,t.power2)
  write.csv(t.power.array2[, , b], file = paste0(dir, paste0("knownnew_t_power_rho", cov_rho, "signalrunif_beta_", betashare[b],".csv")))
}

###############################################
