library(GPCERF)
library(SuperLearner)
library(ggplot2)

set.seed(781)
m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = 12, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

w_all = c(10)
params_lst_50 <- list(alpha = 1, beta = 1, g_sigma = 1,
                      tune_app = "all", n_neighbor = 50, block_size = 1e3)
params_lst_100 <- list(alpha = 1, beta = 1, g_sigma = 1,
                       tune_app = "all", n_neighbor = 100, block_size = 1e3)
params_lst_200 <- list(alpha = 1, beta = 1, g_sigma = 1,
                       tune_app = "all", n_neighbor = 200, block_size = 1e3)

# File to store runtime data
runtime_file <- "runtime_data_var_nn_2.txt"
file.create(runtime_file)

lapply(seq(5000, 100000, 5000), function(n){
  print(n)
  sim_data <- generate_synthetic_data(sample_size = n, gps_spec = 1)
  gps_m <- estimate_gps(cov_mt = sim_data[,-(1:2)],
                        w_all = sim_data$treat,
                        sl_lib = c("m_xgboost", "m_ranger"),
                        dnorm_log = TRUE)

  t1 <- system.time(cerf_nngp_obj <- estimate_cerf_nngp(
                                       sim_data,
                                       w_all,
                                       gps_m,
                                       params = params_lst_50,
                                       outcome_col = "Y",
                                       treatment_col = "treat",
                                       covariates_col = paste0("cf", seq(1,6)),
                                       nthread = 1))
  t2 <- system.time(cerf_nngp_obj <- estimate_cerf_nngp(
                                       sim_data,
                                       w_all,
                                       gps_m,
                                       params = params_lst_100,
                                       outcome_col = "Y",
                                       treatment_col = "treat",
                                       covariates_col = paste0("cf", seq(1,6)),
                                       nthread = 1))
  t3 <- system.time(cerf_nngp_obj <- estimate_cerf_nngp(
                                       sim_data,
                                       w_all,
                                       gps_m,
                                       params = params_lst_200,
                                       outcome_col = "Y",
                                       treatment_col = "treat",
                                       covariates_col = paste0("cf", seq(1,6)),
                                       nthread = 1))

  # Write to file
  cat(paste(n, "nnGP", 50, t1[3]), file = runtime_file,
      append = TRUE, sep = "\n")
  cat(paste(n, "nnGP", 100, t2[3]), file = runtime_file,
      append = TRUE, sep = "\n")
  cat(paste(n, "nnGP", 200, t3[3]), file = runtime_file,
      append = TRUE, sep = "\n")
})
