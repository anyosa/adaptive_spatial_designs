# args = commandArgs(trailingOnly=TRUE)

library(parallel)

source('functions.R')
load('objects.RData')

number_of_designs = 4 
print(paste('Number of designs is:', number_of_designs))

coords_grid = grid(n_northing = 5, n_easting = 5, unitsq = FALSE) 
coords_unitsq = grid(n_northing = 5, n_easting = 5, unitsq = TRUE) 

# number of iterations for the mcmc approximation
iter = 1000000

# we ran this function using parallel, for each k prior scenario, defined by the parameters in matrix_params

get_results <- function(k){
  
  priors <- get_priors(matrix_params, k, coords_unitsq)
  mu_prior <- priors$mu_prior
  Sigma_prior <- priors$Sigma_prior
  
  EIBV_star_designs = c()
  time_star = c()
  
  for(d in 1:number_of_designs){
    eval_d = designs[[d]]
    start = proc.time()
    EIBV_star_designs[d] = EIBV_star_test(coords_grid = coords_grid, mu_prior = mu_prior, Sigma_prior = Sigma_prior, eval_d = eval_d)
    end = proc.time()
    time_star[d] <- (end-start)['elapsed']
    
  }
  
  EIBV_dagger_designs = c() 
  time_dagger = c()
  
  for(d in 1:number_of_designs){
    eval_d = designs[[d]]
    start = proc.time()
    EIBV_dagger_designs[d] = EIBV_dagger_test(coords_grid = coords_grid, mu_prior = mu_prior, Sigma_prior = Sigma_prior, B = 1000, iter_gaussian_approx = 15, eval_d = eval_d)$EIBV_dagger
    end = proc.time()
    time_dagger[d] <- (end-start)['elapsed']
  }
  
  EIBV_MCMC_designs = c() 
  time_mcmc = c()
  
  for(d in 1:number_of_designs){
    eval_d = designs[[d]]
    start = proc.time()
    EIBV_MCMC_designs[d] = EIBV_mcmc_test(coords_grid, mu_prior, Sigma_prior, B = 1000, ITERMAX = iter, tuning = 0.015, eval_d = eval_d)$EIBV_MCMC
    end = proc.time()
    time_mcmc[d] <- (end-start)['elapsed']
  }
  
  return(list(eibv_star =  EIBV_star_designs, eibv_dagger = EIBV_dagger_designs, eibv_mcmc = EIBV_MCMC_designs,
              time_star = time_star, time_dagger = time_dagger, time_mcmc = time_mcmc))
}

numCores <- 4 # args[1]
print(paste('number of cores is:', numCores))

# sequence of prior scenarios to run:
k <- seq(1, 100)

start = proc.time()
results <- mclapply(k, get_results, mc.cores = numCores) # results: a list containing the eibv values and times.
end = proc.time()
print(end-start)