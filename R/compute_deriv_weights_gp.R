#' @title
#' Calculate Derivatives of CERF
#'
#' @description
#' Calculate the weights assigned to each observed outcome when deriving the
#' posterior mean of the first derivative of CERF at a given exposure level.
#'
#' @param w A scalar of exposure level of interest.
#' @param w_obs A vector of observed exposure levels of all samples.
#' @param GPS_m A data.table of GPS vectors.
#'   - Column 1: GPS
#'   - Column 2: Prediction of exposure for covariate of each data sample (e_gps_pred).
#'   - Column 3: Standard deviation of  e_gps (e_gps_std)
#' @param hyperparam A vector of hyper-parameters in the GP model.
#' @param kernel_fn The covariance function.
#' @param kernel_deriv_fn The partial derivative of the covariance function.
#'
#' @return
#' A vector of weights for all samples, based on which the posterior mean of the derivative of CERF at the
#' exposure level of interest is calculated.
#' @export
#'
#' @examples
#'
#' set.seed(915)
#' data <- generate_synthetic_data(sample_size = 200)
#' GPS_m <- train_GPS(cov.mt = as.matrix(data[,-(1:2)]),
#'                    w.all = as.matrix(data$treat))
#'
#' wi <- 4.8
#' weights <- compute_deriv_weights_gp(w = wi,
#'                                     w_obs = data$treat,
#'                                     GPS_m = GPS_m,
#'                                     hyperparam = c(1,1,2))
#'
compute_deriv_weights_gp <- function(w,
                                     w_obs,
                                     GPS_m,
                                     hyperparam,
                                     kernel_fn = function(x) exp(-x),
                                     kernel_deriv_fn = function(x) -exp(-x)){


  alpha <- hyperparam[[1]]
  beta <- hyperparam[[2]]
  g_sigma <- hyperparam[[3]]


  GPS <- GPS_m$GPS
  e_gps_pred <- GPS_m$e_gps_pred
  e_gps_std <- GPS_m$e_gps_std


  # param[1]: alpha, param[2]: beta, param[3]: gamma
  # cov = gamma*h(alpha*w^2 + beta*GPS^2) + diag(1)
  GPS_w <- dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = TRUE)
  n <- length(GPS_w)

  obs_use <- cbind( w_obs*sqrt(1/alpha), GPS*sqrt(1/beta) )
  obs_new <- cbind( w*sqrt(1/alpha), GPS_w*sqrt(1/beta) )
  Sigma_obs <- g_sigma*kernel_fn(as.matrix(dist(obs_use))^2) + diag(nrow(obs_use))
  cross_dist <- spatstat.geom::crossdist(obs_new[,1], obs_new[,2],
                                        obs_use[,1], obs_use[,2])
  Sigma_cross <- g_sigma*sqrt(1/alpha)*kernel_deriv_fn(cross_dist^2)*
    (2*outer(rep(w,n), w_obs, "-"))
  weights_all <- Sigma_cross%*%chol2inv(chol(Sigma_obs))
  # weights_all[weights_all<0] = 0
  # weights = colMeans(weights_all)
  # weights/sum(weights)
  return(colMeans(weights_all))
}
