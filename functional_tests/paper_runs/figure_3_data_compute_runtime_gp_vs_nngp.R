library(GPCERF)
library(SuperLearner)
library(ggplot2)

# Function definitions
m_xgboost <- function(nthread = 4, ...) {
  SuperLearner::SL.xgboost(nthread = nthread, ...)
}

m_ranger <- function(num.threads = 4, ...){
  SuperLearner::SL.ranger(num.threads = num.threads, ...)
}

# Parameters
w_all = c(10)
params_lst <- list(alpha = 1,
                   beta = 1,
                   g_sigma = 1,
                   tune_app = "all",
                   n_neighbor = 50,
                   block_size = 1e3)

# File to store runtime data
runtime_file <- "runtime_data_3.txt"
file.create(runtime_file)

# Generate a vector of 20 random seed values
# set.seed(568)
# seed_values <- sample(10000:20000, 20, replace = FALSE)
seed_values <- c(19259, 12216, 17973, 10706, 14989, 18401, 16808,
                 16002, 16605, 12673, 13626, 16085, 13846, 18507,
                 10682, 18505, 14516, 19595, 13994, 17511)

for (seed in seed_values) {
  set.seed(seed)
  lapply(seq(3000, 10000, 1000), function(n){
    print(paste("Processing sample size:", n, "with seed:", seed))
    sim_data <- generate_synthetic_data(sample_size = n, gps_spec = 1)
    gps_m <- estimate_gps(cov_mt = sim_data[,-(1:2)], w_all = sim_data$treat,
                          sl_lib = c("m_xgboost", "m_ranger"), dnorm_log = TRUE)

    t1 <- system.time(estimate_cerf_gp(sim_data, w_all, gps_m,
                                       params = params_lst,
                                       outcome_col = "Y",
                                       treatment_col = "treat",
                                       covariates_col = paste0("cf", seq(1,6)),
                                       nthread = 1))
    t2 <- system.time(estimate_cerf_nngp(sim_data,
                                         w_all,
                                         gps_m,
                                         params = params_lst,
                                         outcome_col = "Y",
                                         treatment_col = "treat",
                                         covariates_col = paste0("cf", seq(1,6)),
                                         nthread = 1))

    # Write to file
    cat(paste(seed, n, "GP", t1[3]), file = runtime_file,
        append = TRUE, sep = "\n")
    cat(paste(seed, n, "nnGP", t2[3]), file = runtime_file,
        append = TRUE, sep = "\n")
  })
}
