test_that("estimate_cerf_nngp works as expected!", {

  set.seed(19)
  sim.data <- generate_synthetic_data(sample_size = 200, gps_spec = 3)
  # Estimate GPS function
  GPS_m <- train_GPS(cov.mt = as.matrix(sim.data[,-(1:2)]),
                     w.all = as.matrix(sim.data$treat))
  # exposure values
  w.all <- seq(0,20,0.5)
  data.table::setDT(sim.data)
  cerf_nngp_obj <- estimate_cerf_nngp(sim.data,
                                      w.all,
                                      GPS_m,
                                      params = list(alpha = c(0.1,0.2),
                                                    beta = 0.2,
                                                    g_sigma = 1,
                                                    tune_app = "all",
                                                    n_neighbor = 20,
                                                    expand = 1,
                                                    block_size = 1e4))

  expect_s3_class(cerf_nngp_obj, "cerf_nngp")

  expect_equal(length(cerf_nngp_obj$pst_mean), 41L)
  expect_equal(length(cerf_nngp_obj$w), 41L)
  expect_equal(cerf_nngp_obj$pst_mean[1], -22.21065, tolerance = 0.00001)
  expect_equal(cerf_nngp_obj$pst_mean[10], -8.955434, tolerance = 0.00001)
  expect_equal(cerf_nngp_obj$w[31], w.all[31], tolerance = 0.00001)

})
