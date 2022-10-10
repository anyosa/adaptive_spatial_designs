# CSDA #

library(mvtnorm)
library(magrittr)
library(dplyr)
library(tidyr)

generate_all_data <- function(eval_d){
  size_d <- eval_d$size
  list_var <- list()
  for(i in 1:size_d){
    list_var[[i]] <- 0:1
  }
  df <- arrange_all(expand.grid(list_var))
  return(df)
}

filter_with_all_data <- function(y_d, all_data){
  y_d <- matrix(y_d)
  if(nrow(y_d) == ncol(all_data)){
    row_number <- which(apply(all_data, MARGIN = 1, function(x) {all(x == y_d)}) == TRUE)
    return(row_number)
  }
  else(return('error'))
}

grid <- function(n_northing, n_easting, unitsq = TRUE){
  coords <- expand.grid(1:n_northing, 1:n_easting)
  names(coords) <- c("northing", "easting")
  if(unitsq == TRUE){
    coords$northing = coords$northing/n_northing
    coords$easting = coords$easting/n_easting
  } else{
    coords$northing = coords$northing
    coords$easting = coords$easting
  }
  return(coords)
}

matern <- function(zeta, phi, coords) {
  h_matrix = as.matrix(dist(coords))
  cov = (zeta**2) * (1 + phi * h_matrix) * exp(-phi * h_matrix)
  cov
}

simulate_normal <- function(mu, Sigma){
  n = nrow(mu)
  simulated_normal = mu + t(chol(Sigma))%*% rnorm(n)
  return(simulated_normal)
}

simulate_y <- function(linear_predictor){
  m = nrow(linear_predictor)
  eta = linear_predictor
  prob = plogis(eta)
  y = vector("numeric")
  for(i in 1:m){y[i] = rbinom(1,1,prob[i])}
  return(as.matrix(y, nrow = m))
}

selection_matrix <- function(n, index_d){
  size_d = length(index_d)
  A = matrix(data = 0L, nrow = size_d, ncol = n)
  for(i in 1:size_d){
    A[i,index_d[i]] = 1
  }
  return(A)
}

index_design <- function(coords, design){
  size_d = design$size 
  D = data.frame()
  index_d = c()
  for(i in 1:size_d){
    D = rbind(D, filter(coords, northing %in% design$coords_northing[i] & easting %in% design$coords_easting[i]))
    index_d = append(index_d, which(coords$northing == design$coords_northing[i] & coords$easting == design$coords_easting[i]))
  }
  return(index_d)
}


kappa2 <- function(eta_star){
  value = ((1+exp(eta_star))**2) / exp(eta_star)
  return(value)
}

transformed_z <- function(y, eta_star){
  sq_part = (1 + exp(eta_star))**2
  value = ( y - plogis(eta_star) + eta_star * exp(eta_star) / sq_part ) / (exp(eta_star) / sq_part)
  return(value)
}

mu_hat <- function(mu, Sigma_point, Sigma_d, K_d, z_d, mu_d){
  value = mu + Sigma_point %*% solve(Sigma_d + K_d) %*% (z_d - mu_d)
  return(value)
}

Sigma_hat <- function(Sigma, Sigma_point, Sigma_d, K_d){
  value = Sigma - Sigma_point %*% solve(Sigma_d +K_d) %*% t(Sigma_point)
  return(value)
}

chisq <- function(Sigma_point, Sigma_d, K_d){
  matrix_for_d = Sigma_point %*% solve(Sigma_d + K_d) %*% t(Sigma_point)
  value = diag(matrix_for_d)
  return(value)
}

langevin_sampler = function(n, y_d, A, mu_prior, Sigma_prior, ITERMAX, tuning){
  
  mm = mu_prior 
  SS = Sigma_prior 
  L = chol(SS)
  Q = solve(SS)
  
  x = matrix(rep(0, n*ITERMAX), nrow = n, ncol = ITERMAX) 
  eta_C = mu_prior 
  x[, 1]= eta_C
  eta_C_d = A %*% eta_C
  
  logtargC = sum(eta_C_d*y_d - log(1+exp(eta_C_d)))-0.5*as.numeric(t(eta_C-mm)%*%Q%*%(eta_C-mm))
  gradient_C = complete_y_vector(A, y_d - plogis(eta_C_d)) - Q%*%(eta_C-mm) 
  iter = 1
  yes_rate = 0
  
  while(iter < ITERMAX){
    eta_P = eta_C + tuning*gradient_C + sqrt(2*tuning)*rnorm(n)
    eta_P_d = A %*% eta_P
    logtargP = sum(eta_P_d*y_d - log(1+exp(eta_P_d)))-0.5*as.numeric(t(eta_P-mm)%*%Q%*%(eta_P-mm))
    gradient_P = complete_y_vector(A, y_d - plogis(eta_P_d)) - Q%*%(eta_P-mm)
    
    
    q_P = as.numeric(-0.5*(1/(2*tuning))*t(eta_P - eta_C - tuning*gradient_C)%*%(eta_P-eta_C-tuning*gradient_C))
    q_C = as.numeric(-0.5*(1/(2*tuning))*t(eta_C - eta_P - tuning*gradient_P)%*%(eta_C-eta_P-tuning*gradient_P))
    alf = min(1, exp(logtargP-logtargC + q_C - q_P))
    if(alf > runif(1)){
      x[,iter+1] = eta_P
      eta_C = eta_P
      logtargC=logtargP
      gradient_C = gradient_P
      yes_rate = yes_rate + 1
      
    } else{ 
      x[,iter+1] = x[,iter]}
    iter = iter+1
  }
  rate = yes_rate/ITERMAX
  
  take_each <- function(j){
    each <- 1000
    chain <- x[j, ]
    chain_ <- matrix(chain, ncol = each, byrow = TRUE)[,each]
    return(chain_)
  }
  samples_after_thinning <- t(sapply(seq(1:n), take_each))
  return(list(samples = samples_after_thinning, rate = rate, all_samples = x))
}

EBV <- function(u, q){
  correlation = diag(2)
  correlation[lower.tri(correlation)] = q
  correlation[upper.tri(correlation)] = q
  value = pmvnorm(lower = -Inf, upper = c(u, -u), corr = correlation)
  return(value[1])
}

gaussian_approximation <- function(y_d, A, mu_prior, Sigma_prior, iter){ 
  mu_d = A %*% mu_prior
  Sigma_d = A %*% Sigma_prior %*% t(A)
  Sigma_point = Sigma_prior %*% t(A)
  eta_star_d = mu_d
  kappa2_value_d = kappa2(eta_star_d)
  if(nrow(kappa2_value_d) == 1){K_d = kappa2_value_d } #dim 1 1
  if(nrow(kappa2_value_d) != 1){K_d = diag(as.vector(kappa2_value_d))}
  for(i in 1:iter){
    z_d = transformed_z(y_d, eta_star_d)
    mu_posterior = mu_hat(mu = mu_prior, Sigma_point = Sigma_point, Sigma_d = Sigma_d, K_d = K_d, z_d = z_d, mu_d = mu_d) #ok
    eta_star_d = A %*% mu_posterior
    kappa2_value_d = kappa2(eta_star_d)
    if(nrow(kappa2_value_d) == 1){K_d = kappa2_value_d } 
    if(nrow(kappa2_value_d) != 1){K_d = diag(as.vector(kappa2_value_d))}
  }
  Sigma_posterior = Sigma_hat(Sigma = Sigma_prior, Sigma_point = Sigma_point, Sigma_d = Sigma_d, K_d = K_d)
  return(list(mu = mu_posterior, Sigma = Sigma_posterior, z_d = z_d, K_d = K_d)) # modified May28
}

logistic_approximation <- function(mean, cov){
  p = pnorm((0.58*mean)/(sqrt(1+(0.58**2)*diag(cov))))
  return(p)
}

IBV_out_of_d <- function(probabilities, A){
  phis = probabilities
  EBV_s = phis*(1-phis)
  EBV_d = A %*% EBV_s
  return( sum(EBV_s) - sum(EBV_d) )
}

EIBV_star_test <- function(coords_grid, mu_prior, Sigma_prior, eval_d){
  
  n = nrow(coords_grid)
  index_d = index_design(coords = coords_grid, design = eval_d)
  A = selection_matrix(n = n, index_d) 

  eta_star = mu_prior 
  kappa2_value = kappa2(eta_star)
  K = diag(as.vector(kappa2_value))
  K_d = A %*% K %*% t(A)
  mu_d = A %*% mu_prior
  Sigma_d = A %*% Sigma_prior %*% t(A)
  Sigma_point = Sigma_prior %*% t(A)
  
  Sigma_posterior = Sigma_hat(Sigma = Sigma_prior, Sigma_point = Sigma_point, Sigma_d = Sigma_d, K_d = K_d)
  
  alpha_value = 0.58
  sigma2_posterior = diag(Sigma_posterior)
  Chi2 = chisq(Sigma_point = Sigma_point, Sigma_d = Sigma_d, K_d = K_d)
  # r_values = r(alpha_value = alpha_value, sigma2 = sigma2_posterior)
  r_values = alpha_value / (sqrt(1 + alpha_value**2 * sigma2_posterior))
  u_values = (r_values * mu_prior) / sqrt(1 + r_values**2 * Chi2)
  q_values = (-Chi2 * r_values**2) / (1+r_values**2 * Chi2)
  
  EBV_s = c()
  for(i in 1:n){
    EBV_s[i] = EBV(u_values[i], q_values[i])
  }
  EBV_d = A %*% EBV_s 
  EIBV_star = sum(EBV_s) - sum(EBV_d)
  return(EIBV_star)
}

EIBV_dagger_test <- function(coords_grid, mu_prior, Sigma_prior, B, iter_gaussian_approx, eval_d){
  
  n = nrow(coords_grid)
  all_data = generate_all_data(eval_d)
  index_d = index_design(coords = coords_grid, design = eval_d)
  A = selection_matrix(n = n, index_d) 
  EBV_all_data_d = c()
  for(k in 1:NROW(all_data)){
    y_d = as.numeric(all_data[k,])
    params = gaussian_approximation(y_d = y_d, A, mu_prior, Sigma_prior, iter = iter_gaussian_approx)
    phis = logistic_approximation(mean = params$mu , cov = params$Sigma)
    EBV_all_data_d[k] = IBV_out_of_d(probabilities = phis, A)
  }
  
  ws = c()
  list_y_d <- list()
  for(b in 1:B){
    set.seed(b) 
    eta = simulate_normal(mu = mu_prior, Sigma = Sigma_prior)
    eta_d = A %*% eta
    y_d = simulate_y(linear_predictor = eta_d) # nx1
    list_y_d[[b]] <- y_d
    kth_value = filter_with_all_data(y_d, all_data)
    ws[b] = EBV_all_data_d[kth_value]
  }
  
  EIBV_dagger = mean(ws)
  return(list(EIBV_dagger = EIBV_dagger, EBV_all_data_d = EBV_all_data_d, list_y_d = list_y_d))
}

EIBV_mcmc_test <- function(coords_grid, mu_prior, Sigma_prior, B, ITERMAX, tuning = 0.0018, eval_d){
  
  n = nrow(coords_grid)
  all_data = generate_all_data(eval_d)
  index_d = index_design(coords = coords_grid, design = eval_d)
  A = selection_matrix(n = n, index_d) 
  matrix_y_d = matrix(rep(0,eval_d$size*B), ncol = eval_d$size)  
  list_y_d <- list()
  
  for(b in 1:B){
    set.seed(b)
    eta = simulate_normal(mu_prior, Sigma_prior)
    eta_d = A %*% eta
    y_d = simulate_y(linear_predictor = eta_d)
    matrix_y_d[b,] = y_d
    list_y_d[[b]] <- y_d
    df = as.data.frame(matrix_y_d)
  }
  combinations = distinct(df)
  
  EBV_in_combinations = c()
  for(r in 1:nrow(combinations)){
    y_d = as.numeric(combinations[r,])
    samples = langevin_sampler(n, y_d, A, mu_prior, Sigma_prior, ITERMAX, tuning)$sample # n x thinned_samples
    EBV_in_combinations[r] = EBV_MCMC(samples_matrix = samples, A)
  }
  
  temp = c()
  ws = c()
  for(b in 1:B){
    y_d = matrix_y_d[b,]
    kth_value = filter_with_all_data(y_d, combinations)
    temp[b] = kth_value
    ws[b] = EBV_in_combinations[kth_value]
  }
  EIBV_MCMC = mean(ws)
  return(list(EIBV_MCMC = EIBV_MCMC, list_y_d = list_y_d))
}

complete_y_vector <- function(A, vector_d){
  index_in_y <- which(A==1, arr.ind = TRUE)[,2]
  y_zeros = rep(0, ncol(A))
  size_d = nrow(A)
  for(i in 1:size_d){
    y_zeros[index_in_y[i]] = vector_d[i]
  }
  recover <- A%*%y_zeros
  if(sum(recover == vector_d) == size_d){
    return(y_zeros)
  } else{print('no match')}
}

EBV_MCMC = function(samples_matrix, A){
  probs_matrix = plogis(samples_matrix) 
  probs_hat = rowMeans(probs_matrix) 
  EBV_s = probs_hat*(1-probs_hat)
  EBV_d = A %*% EBV_s 
  return(sum(EBV_s) - sum(EBV_d))
}

set_design = function(coords_northing, coords_easting, n_northing, n_easting){
  if(length(coords_northing)==length(coords_easting)){
    output = list(
      coords_northing = coords_northing,
      coords_easting = coords_easting,
      n1 = n_northing,
      n2 = n_easting,
      size = length(coords_northing))
  }else{print('Error')}
}

generate_features <- function(coords, northing_center, easting_center){
  feature <- function(coords){(coords$northing-northing_center)**2 + (coords$easting-easting_center)**2}
  X <- matrix(cbind(rep(1,NROW(coords)), feature(coords)), nrow = nrow(coords))
  return(X)
}

get_priors <- function(matrix_params, k, coords_unitsq){
  vec_1 <- matrix_params[k, 1]
  vec_2 <- matrix_params[k, 2]
  X <- generate_features(coords_unitsq, vec_1/10, vec_2/10) 
  n = nrow(X)
  
  mu_1 <- matrix_params[k, 3]
  mu_2 <- matrix_params[k, 4]
  mu_beta = matrix(c(mu_1, mu_2), 2, 1)
  
  sd_1 <- matrix_params[k, 5]
  sd_2 <- matrix_params[k, 6]
  rho <- matrix_params[k, 7]
  Sigma_beta = matrix(c(sd_1**2, rho*sd_1*sd_2, rho*sd_1*sd_2, sd_2**2), 2, 2)
  
  mu_w = matrix(rep(0, n), n, 1)
  zeta <- matrix_params[k, 8]
  phi <- matrix_params[k, 9]
  Sigma_w = matern(zeta = zeta, phi = phi, coords = coords_unitsq)
  
  params <- c(vec_1, vec_2, mu_1, mu_2, sd_1, sd_2, rho, zeta, phi)
  
  mu_prior = X %*% mu_beta
  Sigma_prior = X %*% Sigma_beta %*% t(X) + Sigma_w
  return(list(mu_prior = mu_prior, Sigma_prior = Sigma_prior))
}

get_possible_points = function(matrix_current_set_points, number_rows_grid, number_cols_grid){
  point = tail(matrix_current_set_points, 1)
  possible = as.data.frame(
    matrix(c(point[1] - 1, point[2],
             point[1] + 1, point[2],
             point[1], point[2] - 1,
             point[1], point[2] + 1), byrow = TRUE, ncol = 2))
  
  possible = filter(possible, !(V1 %in% c(0, number_rows_grid+1) | V2 %in% c(0, number_cols_grid+1) ))
  k = nrow(matrix_current_set_points)
  if(k>1){
    N = k-1
    previous = head(matrix_current_set_points, N)
    for(n in 1:N){
      possible = filter(possible, !(V1 %in% previous[n,1] & V2 %in% previous[n,2])) # remove previous point from possible ones
    }
  }
  return(as.matrix(possible))
}

sequential_sampling_AUV_deterministic =  function(number_of_designs, size_of_design, initial_point, mu_prior, Sigma_prior, y, coords_grid){
  
  n = nrow(coords_grid)
  n_northing = max(coords_grid$northing)
  n_easting = max(coords_grid$easting)
  B = number_of_designs
  W = size_of_design - 1
  designs_over_B = list()
  EIBV_over_W_B = matrix(0L, nrow = B, ncol = W+1)
  matrix_IBV = matrix(0L, nrow = B, ncol = W+1)
  
  for(b in 1:B){ 
    eval_d = set_design(initial_point[1], initial_point[2], n_northing, n_easting)
    index_d = index_design(coords = coords_grid, design = eval_d)
    A = selection_matrix(n = n, index_d)
    y_d = A %*% y
    updated_parameters = gaussian_approximation(y_d, A, mu_prior, Sigma_prior, iter = 15)
    mu_prior = updated_parameters$mu
    Sigma_prior =  updated_parameters$Sigma 
    matrix_IBV[b,1] = IBV_out_of_d(probabilities = logistic_approximation(mu_prior, Sigma_prior), A)
    for(w in 1:W){ 
      EIBV_star_points = c()
      if(w == 1){
        possible_points = get_possible_points(initial_point, n_northing, n_easting)
        current_set_points = initial_point
      }
      P = nrow(possible_points) 
      for(p in 1:P){ 
        matrix_points = rbind(current_set_points, possible_points[p,])
        eval_d = set_design(matrix_points[,1], matrix_points[,2], n_northing, n_easting) 
        EIBV_star_points[p] = EIBV_star_test(coords_grid = coords_grid, mu_prior = mu_prior, Sigma_prior = Sigma_prior, eval_d = eval_d)
      }
      index_of_selected_point = match(min(EIBV_star_points), EIBV_star_points) 
      eval_d = set_design(as.numeric(possible_points[index_of_selected_point,][1]), as.numeric(possible_points[index_of_selected_point,][2]), n_northing, n_easting)
      index_d = index_design(coords = coords_grid, design = eval_d)
      A = selection_matrix(n = n, index_d)
      y_d = A %*% y
      updated_parameters = gaussian_approximation(y_d, A, mu_prior, Sigma_prior, iter = 15) 
      mu_prior = updated_parameters$mu 
      Sigma_prior =  updated_parameters$Sigma 
      current_set_points = rbind(current_set_points, possible_points[index_of_selected_point,])
      eval_d = set_design(current_set_points[,1], current_set_points[,2], n_northing, n_easting) 
      index_d = index_design(coords = coords_grid, design = eval_d)
      A = selection_matrix(n = n, index_d)
      matrix_IBV[b,w+1] = IBV_out_of_d(probabilities = logistic_approximation(mu_prior, Sigma_prior), A)
      possible_points = get_possible_points(current_set_points, n_northing, n_easting)
      if(dim(possible_points)[1] == 0){break}
    }
    designs_over_B[[b]] = current_set_points
  }
  return(list(designs = designs_over_B, IBV = matrix_IBV))
}

evaluate_prescripted_designs <- function(design, mu_prior, Sigma_prior, coords_grid, y){
  eval_d <- design
  index_d <- index_design(coords = coords_grid, design = eval_d)
  A <- selection_matrix(n = n, index_d)
  y_d <- A %*% y
  updated_parameters = gaussian_approximation(y_d, A, mu_prior, Sigma_prior, iter = 15)
  mu_post = updated_parameters$mu
  Sigma_post =  updated_parameters$Sigma 
  value_ <- IBV_out_of_d(probabilities = logistic_approximation(mu_post, Sigma_post), A)
  return(value_)
}

evaluate_adaptive_desings <- function(matrix_one_adaptive_path, mu_prior, Sigma_prior, coords_grid, y){
  n_northing <- max(coords_grid$northing)
  n_easting <- max(coords_grid$easting)
  matrix_design <- matrix_one_adaptive_path
  eval_d <- set_design(matrix_design[, 1], matrix_design[, 2], n_northing, n_easting)
  index_d <- index_design(coords = coords_grid, design = eval_d)
  A <- selection_matrix(n = n, index_d)
  y_d <- A %*% y
  updated_parameters = gaussian_approximation(y_d, A, mu_prior, Sigma_prior, iter = 15)
  mu_post = updated_parameters$mu 
  Sigma_post =  updated_parameters$Sigma 
  value_ <- IBV_out_of_d(probabilities = logistic_approximation(mu_post, Sigma_post), A)
  return(value_)
}