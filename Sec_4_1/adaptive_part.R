source('functions.R')
load('objects.RData')

#### (1) Set coordinates. ####

coords_grid = grid(n_northing = 5, n_easting = 5, unitsq = FALSE) 
coords_unitsq = grid(n_northing = 5, n_easting = 5, unitsq = TRUE) 

n = nrow(coords_grid)

#### (2) Work with pre-scripted designs. ####

doe <- list()
doe[[1]] <- set_design(seq(1,5), rep(5,5), 5, 5) 
doe[[2]] <- set_design(rep(5,5), seq(1,5), 5, 5) 

#### (3) Generate adaptive paths. ####

# set1_matrix_ys is the first set of ground truths, we use this in the comparison to d=5
# set2_matrix_ys is the second set of ground truths, we use this in the comparison to d=6
# set1_matrix_initial_points is a matrix with initial points randomly chosen, we use this in the comparison to d=5
# set2_matrix_initial_points is a matrix with initial points randomly chosen, we use this in the comparison to d=6

matrix_ys <- set1_matrix_ys # or set2_matrix_ys
matrix_initial_points <- set1_matrix_initial_points # or set2_matrix_initial_points
start <- proc.time()
K <- nrow(matrix_params)

list_initial_points <- list()
AUV_designs <- list()
AUV_IBV <- list()
for(k in 1:K){
  print(paste('Generating path for replicate:', k))
  priors <- get_priors(matrix_params, k, coords_unitsq)
  mu_prior <- priors$mu_prior
  Sigma_prior <- priors$Sigma_prior
  initial_point <- matrix(matrix_initial_points[k, ], ncol = 2) 
  list_initial_points[[k]] <- initial_point
  y <- matrix(matrix_ys[k, ], nrow = 25)
  AUV_sampling = sequential_sampling_AUV_deterministic(number_of_designs = 1, size_of_design = 5, initial_point, mu_prior, Sigma_prior, y, coords_grid)
  AUV_IBV[[k]] = AUV_sampling$IBV[1, ]
  AUV_designs[[k]] <- AUV_sampling$designs[[1]]
}
end <- proc.time()
(end-start)['elapsed'] 

size_of_design = 5
get_vector_of_IBV <- function(){
  k <- seq(1, 100)
  IBV <- sapply(k, function(k) AUV_IBV[[k]][size_of_design])
  to_fix <- which(sapply(k, function(k) AUV_IBV[[k]][size_of_design]) == 0)
  if(length(to_fix) > 0){
    for(j in 1:length(to_fix)){
      m <- min(which(AUV_IBV[[to_fix[j]]] == 0))-1
      IBV[to_fix[j]] <- AUV_IBV[[to_fix[j]]][m]
    }
    which(IBV == 0)
  }
  return(IBV)
}

IBV_a <- get_vector_of_IBV() # adaptive IBV

#### (4) Compare paths ####

p = 1 # p = 1 is d = 5 and p = 2 is d = 6
start <- proc.time()
IBV_p <- c()
for(k in 1:K){
priors <- get_priors(matrix_params, k, coords_unitsq)
mu_prior <- priors$mu_prior
Sigma_prior <- priors$Sigma_prior
y <- matrix(matrix_ys[k, ], nrow = 25)
IBV_p[k] <- evaluate_prescripted_designs(design = doe[[p]], mu_prior, Sigma_prior, coords_grid, y)
}
end <- proc.time()
(end-start)['elapsed']

# to compare we can use the vectors from results.RData
IBV_a <- set1_IBV_adaptive # or set_2_IBV_adaptive
IBV_p <- set1_IBV_prescripted # or set_2_IBV_prescripted

t.test(IBV_a, IBV_p, paired = T) 
diff <- IBV_a - IBV_p
mean(diff)
sd(diff)/sqrt(100)
mean(IBV_a)
sd(IBV_a)/sqrt(100)
mean(IBV_p) 
sd(IBV_p)/sqrt(100) 
