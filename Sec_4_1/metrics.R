# CSDA

# values
k = 1:100

eibv_star <- t(sapply(k, function(k) results[[k]]$eibv_star))
eibv_dagger <- t(sapply(k, function(k) results[[k]]$eibv_dagger))
eibv_mcmc <- t(sapply(k, function(k) results[[k]]$eibv_mcmc))

time_star <-  t(sapply(k, function(k) results[[k]]$time_star))
time_dagger <-  t(sapply(k, function(k) results[[k]]$time_dagger))
time_mcmc <-  t(sapply(k, function(k) results[[k]]$time_mcmc))

#### spearman ####

rank_correlation <- function(matrix_1, matrix_2, cor_method){
  n <- NROW(matrix_1)
  estimates <- c()
  for(i in 1:n){
    estimates[i] <- cor.test(order(matrix_1[i,]), order(matrix_2[i,]), method = cor_method)$estimate
  }
  mean <- round(mean(estimates), 2)
  se <- round(sd(estimates)/sqrt(n), 2)
  return(list(mean = mean, se = se))
}

rank_correlation(eibv_mcmc, eibv_dagger, 'spearman') 
rank_correlation(eibv_mcmc, eibv_star, 'spearman')
rank_correlation(eibv_dagger, eibv_star, 'spearman')

#### rankings ####

compare_order <- function(matrix_1, matrix_2){
  n <- NROW(matrix_1)
  reps <- NROW(matrix_1)
  vals_1 <- c()
  for(i in 1:reps){
    vals_1[i] <- ifelse(order(matrix_1[i, ])[1] == order(matrix_2[i, ])[1], 1, 0)
  }
  prop <- sum(vals_1)/n
  approx_se <- round(sqrt(prop*(1-prop)/n), 2)
  return(list(p = prop, se = approx_se))
}

compare_order(eibv_mcmc, eibv_dagger)
compare_order(eibv_mcmc, eibv_star)
compare_order(eibv_dagger, eibv_star)

#### bias ####

calculate_bias <- function(matrix_1, matrix_2){
  n <- NROW(matrix_1)
  bias <- as.vector(matrix_2 - matrix_1)
  mean <- round(mean(bias), 2)
  se <- round(sd(bias)/sqrt(n), 2)
  return(list(mean = mean, se = se))
}

calculate_bias(eibv_mcmc, eibv_dagger)
calculate_bias(eibv_mcmc, eibv_star)
calculate_bias(eibv_dagger, eibv_star)

#### RMSE ####

calculate_rmse <- function(data, indices){
  d <- data[indices,]
  diff2 <- (d[,1] - d[,2])**2
  total <- length(diff2)
  result <- round(sqrt(sum(diff2)/total), 2)
  return(result)
}

library(boot)
set.seed(123)
matrix_data <- cbind(as.vector(eibv_mcmc),as.vector(eibv_dagger))
boot(matrix_data, calculate_rmse, R = 100)
matrix_data <- cbind(as.vector(eibv_mcmc),as.vector(eibv_star))
boot(matrix_data, calculate_rmse, R = 100)
matrix_data <- cbind(as.vector(eibv_dagger),as.vector(eibv_star))
boot(matrix_data, calculate_rmse, R = 100)


#### time ####

report_time <- function(matrix_time){
  vec_time <- as.vector(matrix_time)
  n <- length(vec_time)
  mean_global <- round(mean(vec_time), 2)
  se_global <- round(sd(vec_time)/sqrt(n), 2)
  return(list(mean = mean_global, se = se_global))
}

report_time(time_dagger)
report_time(time_star)
report_time(time_mcmc)
