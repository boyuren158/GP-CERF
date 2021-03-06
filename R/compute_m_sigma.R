#' @title
#' Compute mean, confidence interval (?), and covariate balance in Full Gaussian
#'  Process (GP)
#'
#' @description
#' Calculates the induced covariate balance associated with one hyper-parameter
#' configuration in full GP.
#'
#' @param hyperparam A vector of values of hyper-parameters.
#'   - First element: alpha
#'   - Second element: beta
#'   - Third element: g_sigma (gamma/sigma)
#' @param data A  data.table containing all data including outcome, exposure
#' and covariates. In the following order:
#'   - Column 1: Outcome (Y)
#'   - Column 2: Exposure or treatment (w)
#'   - Column 3~m: Confounders (C)
#' @param w A vector of exposure levels at which the CERF is estimated.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: A vector of estimated GPS evaluated at the observed exposure levels.
#'   - Column 2: Estimated conditional means of the exposure given covariates
#'               for all samples (e_gps_pred).
#'   - Column 3: Estimated conditional standard deviation of the exposure given
#'               covariates for all samples (e_gps_std).
#' @param nthread An integer value that represents the number of threads to be
#' used by internal packages.
#' @param kernel_fn The covariance function of GP.
#'
#' @return
#' A list containing two elements: 1) a vector of absolute weighted correlation of each
#' covariate to the exposure, which is the metric for covariate balance and 2) the estimated
#' CERF at \code{w.all} based on the hyper-parameter values in \code{param}.
#' @export
#'
#' @examples
#'
#' set.seed(912)
#' data <- generate_synthetic_data(sample_size = 250, gps_spec = 3)
#'
#' w.all = seq(0,20,1)
#'
#' data.table::setDT(data)
#'
#' #Estimate GPS function
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
#'
#' tune_res <- compute_m_sigma(hyperparam = c(0.09, 0.09, 10),
#'                             data = data,
#'                             w = w.all,
#'                             GPS_m = GPS_m,
#'                             nthread = 1)
#'
#' gp.cerf <- tune_res$est
#'
compute_m_sigma <- function(hyperparam, data, w, GPS_m, nthread = 1,
                            kernel_fn = function(x) exp(-x^2)){

  param = unlist(hyperparam)

  GPS <- GPS_m$GPS
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std

  # mi(w)
  # param 1: alpha
  # param 2: beta
  # param 3: ratio gamma/sigma

  alpha <- param[1]
  beta  <- param[2]
  g_sigma <- param[3]


  w_obs <- data[[2]]

  scaled_obs <- cbind(w_obs*sqrt(1/alpha), GPS*sqrt(1/beta))
  sigma_obs <- g_sigma*kernel_fn(as.matrix(dist(scaled_obs))) + diag(nrow(scaled_obs))

  inv_sigma_obs <- compute_inverse(sigma_obs)

  # Estimate noise
  noise_est <- estimate_noise_gp(hyperparam = hyperparam,
                                 data = data,
                                 GPS = GPS)


  lfp <- get_options("logger_file_path")

  # make a cluster
  cl <- parallel::makeCluster(nthread, type="PSOCK",
                              outfile= lfp)

  # export variables and functions to cluster cores
  parallel::clusterExport(cl=cl,
                          varlist = c("w", "w_obs", "scaled_obs",
                                      "hyperparam", "inv_sigma_obs", "GPS_m",
                                      "scale", "nthread", "noise_est",
                                      "compute_sd_gp", "compute_weight_gp",
                                      "compute_w_corr"),
                          envir=environment())

  col_all_list <- parallel::parLapply(cl,
                                      w,
                                      function(w_instance){

    # compute weights
    weights_final <- compute_weight_gp(w = w_instance,
                                       w_obs = w_obs,
                                       scaled_obs = scaled_obs,
                                       hyperparam = hyperparam,
                                       inv_sigma_obs = inv_sigma_obs,
                                       GPS_m = GPS_m)

    weights_final[weights_final<0] <- 0
    weights_final <- weights_final/sum(weights_final)
    # weigts.final = invers of paranthesis * kappa
    # est is the same as m in the paper.
    est <- data$Y%*%weights_final


    pst_sd <- compute_sd_gp(w = w_instance,
                            scaled_obs = scaled_obs,
                            hyperparam = hyperparam,
                            sigma = noise_est,
                            GPS_m = GPS_m)

    covariate_balance <- compute_w_corr(data, weights_final)
    c(covariate_balance, est, pst_sd)
  })

  # terminate clusters.
  parallel::stopCluster(cl)

  col_all <- do.call(cbind, col_all_list)

  n_confounders <- nrow(col_all) - 2 # est, pst_sd
  est_index <- nrow(col_all) - 1
  pst_sd <- nrow(col_all)

  list(cb = rowMeans(col_all[1:n_confounders,], na.rm = T),
       est = col_all[est_index,],
       pst = col_all[pst_sd, ])
}
