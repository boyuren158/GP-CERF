test_that("estimate_mean_sd_nn works as expected!", {

  set.seed(276)
  data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)

  # Estimate GPS function
  GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
                     w.all = as.matrix(data$treat))

  # Hyperparameter
  hyperparam <- c(0.1, 0.2, 1)
  n_neighbor <- 15
  expand <- 1
  block_size <- 10000

  # compute noise
  noise <- estimate_noise_nn(hyperparam = hyperparam,
                             w_obs = data$treat,
                             GPS_obs = GPS_m$GPS,
                             y_obs = data$Y,
                             n_neighbor = n_neighbor)

  # compute posterior mean and standard deviation for vector of w.
  w <- seq(0,20,1)
  val <- estimate_mean_sd_nn(hyperparam = hyperparam,
                             sigma2 = noise,
                             w_obs = data$treat,
                             w = w,
                             y_obs = data$Y,
                             GPS_m = GPS_m,
                             n_neighbor = n_neighbor,
                             expand = expand,
                             block_size = block_size)

  expect_equal(length(val), 2L)
  # There are length(w) elements in the first item and the second item.
  expect_equal(length(val[[1]]), 21L)
  expect_equal(length(val[[2]]), 21L)

  # First element returns number of neighbors + mean value
  expect_equal(nrow(val[[1]][[10]]), 31L)

  expect_equal(val[[1]][[12]][29,2], 0.0003256218, tolerance = 0.0001)

})
