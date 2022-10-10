# This is an example file to introduce the main functions for section 4.1.
# In the pre-scripted part you can run the approximation functions.
# In the adaptive part you can run the function to generate an adaptive path.

# please, set working directory

# load functions and objects (a list with designs and a matrix of hyperparameters - used in section 4.1)
source('functions.R')
load('objects.RData')

# set coordinates for the grid of size 5x5

coords_grid = grid(n_northing = 5, n_easting = 5, unitsq = FALSE) 
coords_unitsq = grid(n_northing = 5, n_easting = 5, unitsq = TRUE) 

# number of iterations for the mcmc approximation
iter = 10000

#### PRE-SCRIPTED PART ####

# choose k: a set of hyper-parameters to define the priors of eta, k =1,...,100
k = 5
priors <- get_priors(matrix_params, k, coords_unitsq)
mu_prior <- priors$mu_prior
Sigma_prior <- priors$Sigma_prior

# choose d: a design d = 1,...,6
d = 1
eval_d = designs[[d]]

# EIBV I
start = proc.time()
EIBV_dagger =  EIBV_dagger_test(coords_grid = coords_grid, mu_prior = mu_prior, Sigma_prior = Sigma_prior, B = 1000, iter_gaussian_approx = 15, eval_d = eval_d)$EIBV_dagger
end = proc.time()
end-start
EIBV_dagger

# EIBV II
start = proc.time()
EIBV_star = EIBV_star_test(coords_grid = coords_grid, mu_prior = mu_prior, Sigma_prior = Sigma_prior, eval_d = eval_d)
end = proc.time()
end-start
EIBV_star

# EIBV MCMC
start = proc.time()
EIBV_mcmc = EIBV_mcmc_test(coords_grid, mu_prior, Sigma_prior, B = 1000, ITERMAX = iter, tuning = 0.015, eval_d = eval_d)$EIBV_MCMC
end = proc.time()
end-start
EIBV_mcmc

#### ADAPTIVE PART ####

# load scenario
k = 1
priors <- get_priors(matrix_params, k, coords_unitsq)
mu_prior <- priors$mu_prior
Sigma_prior <- priors$Sigma_prior
matrix_ys <- set1_matrix_ys

# generate adaptive path
initial_point <- matrix(sample.int(n = 5, size = 2, replace = TRUE, prob = rep(0.2, 5)), ncol = 2)
y <- matrix(matrix_ys[k, ], nrow = 25)
AUV_sampling <- sequential_sampling_AUV_deterministic(number_of_designs = 1, size_of_design = 5, initial_point, mu_prior, Sigma_prior, y, coords_grid)
AUV_designs <- AUV_sampling$designs[[1]]

# see generated path:
# V1 is northing and V2 is easting.
AUV_designs
