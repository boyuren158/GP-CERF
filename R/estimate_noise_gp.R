#' @title
#' Estimate the Standard Deviation of the Nugget Term in Full Gaussian Process
#'
#' @description
#' Estimates the standard deviations of the nugget term in full GP by
#' calculating the standard deviations of the residuals.
#'
#' @param hyperparam A vector of hyper-parameter values for the full GP.
#' @param sigma_obs The covariance matrix between the observed data points.
#' @param inv_sigma_obs The inverse of the covariance matrix between the observed data points.
#'
#' @return
#' A scalar of estimated standard deviation of the nugget term in full GP.
#'
#' @export
#'
#' @examples
#'
#' set.seed(109)
#' data <- generate_synthetic_data(sample_size = 100, gps_spec = 3)
#'
#' # Estimate GPS function
#' GPS_m <- train_GPS(cov_mt = data[,-(1:2)], w_all = data$treat)
#'
#' hyperparam <- c(0.1, 0.2, 1)
#'
#' noise_est <- estimate_noise_gp(hyperparam, data, GPS_m$GPS)
#'
#'
estimate_noise_gp <- function(data, sigma_obs, inv_sigma_obs){
  # Need to update help docs
  noise <- sd(data$Y - arma_mm(sigma_obs - diag(nrow(sigma_obs)),
                               arma_mm(inv_sigma_obs, data$Y)))

  return(noise)
}
