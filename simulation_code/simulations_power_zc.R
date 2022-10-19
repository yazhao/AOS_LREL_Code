###############################################
# This script simulates both null and 
# alternative hypothesis data and checks the 
# power (which in the null case is the type I
# error). We use a beta defined as a delta *
# [1,...,1,0,...,0], starting with 10% 1's
#
# For convenience we've set sigma^2 = 1
#
# This script is for the Zhong and Chen cases
# for normal and t-distribution errors
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
#betashare = c(0.25, 0.5, 0.75, 1)
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
# Normal Error, Full Signal
###############################################
# storage for all different delta values
normal.power.array = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
normal.power = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
normal.power[,1] = delp
# storing for each of the delta runs
normal.run = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results for the zc comparison normal null
  normal_null_zc <- read.csv(paste0(dir,"normal_null_zc_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    delta = delp * mdf[i]
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "full", epsilon = "normal")
        x.mat = results[[2]]
        y.vec = results[[3]]
        
        # calculating the T statistics
        Tstats = zhongchen2011(X = x.mat, y = y.vec, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
        
        # store in matrix
        normal.run[k,j] = Tstats[2]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      normal.power[j, i+1] = sum((abs(normal.run[,j]) > 1.96)*1)/iter
    }
  }
  normal.power.array[, , b] = nullpower_append(normal_null_zc,normal.power)
  write.csv(normal.power.array[, , b], file = paste0(dir,paste0("normal_power_rho", cov_rho, "signalfullzc_beta_",betashare[b],".csv")))
}
###############################################

###############################################
# Normal Error, Runif Signal
###############################################
# storage for all different delta values
normal.power2.array = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
normal.power2 = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
normal.power2[,1] = delp
# storing for each of the delta runs
normal.run2 = matrix(data = 0, nrow = iter, ncol = length(delp))


for(b in 1:length(betashare)){
  # calling in the results for the zc comparison normal null
  normal_null_zc <- read.csv(paste0(dir,"normal_null_zc_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    delta = delp * mdr[i]
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "runif", epsilon = "normal")
        x.mat = results[[2]]
        y.vec = results[[3]]
        
        # calculating the T statistics
        Tstats = zhongchen2011(X = x.mat, y = y.vec, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
        
        # store in matrix
        normal.run2[k,j] = Tstats[2]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      normal.power2[j, i+1] = sum((abs(normal.run2[,j]) > 1.96)*1)/iter
    }
  }
  normal.power2.array[, , b] = nullpower_append(normal_null_zc,normal.power2)
  write.csv(normal.power2.array[, , b], file = paste0(dir,paste0("normal_power_rho", cov_rho, "signalrunifzc_beta_",betashare[b],".csv")))
}
###############################################

###############################################
# t-distribution error, Full Signal
###############################################
# storage for all different delta values
t.power.array = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
t.power = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
t.power[,1] = delp
# storing for each of the delta runs
t.run = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results for the zc comparison normal null
  t_null_zc <- read.csv(paste0(dir,"t_null_zc_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    delta = delp * mdf[i]
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "full", epsilon = "t")
        x.mat = results[[2]]
        y.vec = results[[3]]
        
        # calculating the T statistics
        Tstats = zhongchen2011(X = x.mat, y = y.vec, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
        
        # store in matrix
        t.run[k,j] = Tstats[2]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      t.power[j, i+1] = sum((abs(t.run[,j]) > 1.96)*1)/iter
    }
  }
  t.power.array[, , b] = nullpower_append(t_null_zc,t.power)
  write.csv(t.power.array[, , b], file = paste0(dir,paste0("t_power_rho", cov_rho, "signalfullzc_beta_",betashare[b],".csv")))
}
###############################################

###############################################
# t-distribution error, Runif Signal
###############################################
# storage for all different delta values
t.power2.array = array(data = NA, dim = c(length(delp) + 1, length(n) + 1, length(betashare)))
t.power2 = matrix(data = NA, nrow = length(delp), ncol = length(n) + 1)
t.power2[,1] = delp
# storing for each of the delta runs
t.run2 = matrix(data = 0, nrow = iter, ncol = length(delp))

for(b in 1:length(betashare)){
  # calling in the results for the zc comparison normal null
  t_null_zc <- read.csv(paste0(dir,"t_null_zc_rho", cov_rho, ".csv"))
  
  for(i in 1:length(n)){
    delta = delp * mdr[i]
    for(j in 1:length(delp)){
      for(k in 1:iter){
        # Setting the seed 
        set.seed(k)
        
        # generating the simulated data
        results = simulate_data(n = n[i], c = c, rho = cov_rho, beta = betashare[b] * n[i], delta = delta[j], signal = "runif", epsilon = "t")
        x.mat = results[[2]]
        y.vec = results[[3]]
        
        # calculating the T statistics
        Tstats = zhongchen2011(X = x.mat, y = y.vec, beta = 0, delta = 0, Sigma = NULL, small.sig = NULL, T0 = FALSE)
        
        # store in matrix
        t.run2[k,j] = Tstats[2]
      }
      # calculating power, which is the proportion of test statistics >1.96 (two-sided)
      t.power2[j, i+1] = sum((abs(t.run2[,j]) > 1.96)*1)/iter
    }
  }
  t.power2.array[, , b] = nullpower_append(t_null_zc,t.power2)
  write.csv(t.power2.array[, , b], file = paste0(dir,paste0("t_power_rho", cov_rho, "signalrunifzc_beta_",betashare[b],".csv")))
}
###############################################