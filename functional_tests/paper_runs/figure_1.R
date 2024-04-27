library(GPCERF)

set.seed(781)
# Generate synthetic data with 500 data samples.
sim_data <- generate_synthetic_data(sample_size = 500,
                                    gps_spec = 1)

# SuperLearner internal libraries' wrapper. (Optional)

m_xgboost <- function(nthread = 12, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}
m_ranger <- function(num.threads = 12, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

# Estimate GPS function
gps_m <- estimate_gps(cov_mt = sim_data[, paste0("cf", seq(1,6))],
                      w_all = sim_data$treat,
                      sl_lib = c("m_xgboost", "m_ranger"),
                      dnorm_log = TRUE)

# exposure values of interest
# We trim the exposure level to satisfy positivity assumption to avoid including
# extreme exposure values.
q1 <- stats::quantile(sim_data$treat, 0.05)
q2 <- stats::quantile(sim_data$treat, 0.95)
w_all <- seq(q1, q2, 1)

# Hyperparameters' range for grid search to find optimal hyperparameters
params_lst <- list(alpha = 10 ^ seq(-2, 2, length.out = 10),
                   beta = 10 ^ seq(-2, 2, length.out = 10),
                   g_sigma = c(0.1, 1, 10),
                   tune_app = "all")

# Estimate exposure response function
cerf_gp_obj <- estimate_cerf_gp(sim_data,
                                w_all,
                                gps_m,
                                params = params_lst,
                                outcome_col = "Y",
                                treatment_col = "treat",
                                covariates_col = paste0("cf", seq(1,6)),
                                nthread = 12)


summary(cerf_gp_obj)
plot(cerf_gp_obj)
