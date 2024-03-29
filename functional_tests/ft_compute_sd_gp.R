 set.seed(284)
 #Generate synthetic data
 data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
 w_obs <- obs_exposure <- data$treat

 # Choose an exposure level to compute CERF
 w = 1.2

 # Define kernel function
 kernel_fn <- function(x) exp(-x^2)

 # Estimate GPS function
 GPS_m <- train_GPS(cov_mt = as.matrix(data[,-(1:2)]),
                    w_all = as.matrix(data$treat))

 GPS <- GPS_m$GPS

 # set hyperparameters
 hyperparam <- c(0.1, 0.4, 1)
 alpha <- hyperparam[[1]]
 beta <- hyperparam[[2]]
 g_sigma <- hyperparam[[3]]

 # Compute scaled observation data and inverse of covariate matrix.
 scaled_obs <- cbind(obs_exposure*sqrt(1/alpha), GPS*sqrt(1/beta))

 tentative_sigma <- 0.1

 post_sd <- GPCERF:::compute_sd_gp(w = w,
                                  scaled_obs = scaled_obs,
                                  hyperparam = hyperparam,
                                  sigma = tentative_sigma,
                                  GPS_m = GPS_m,
                                  kernel_fn = kernel_fn)
