#' @title
#' Estimate a model for generalized propensity score
#'
#' @description
#' Estimates a model for generalized propensity score (GPS) using parametric
#' approach.
#'
#' @param cov_mt A covariate matrix containing all covariates. Each row is a
#' data sample and each column is a covariate.
#' @param w_all A vector of observed exposure levels.
#' @param sl_lib A vector of SuperLearner's package libraries.
#' @param dnorm_log Logical, if TRUE, probabilities p are given as log(p).
#'
#' @return
#' A data.frame that includes:
#'   - a vector of estimated GPS at the observed exposure levels;
#'   - a vector of estimated conditional means of exposure levels when the covariates are fixed
#' at the observed values;
#'   - estimated standard deviation of exposure levels
#'   - a vector of observed exposure levels.
#' @export
#'
#' @examples
#' \donttest{
#' data <- generate_synthetic_data(sample_size = 200)
#' GPS_m <- estimate_gps(cov_mt = data[,-(1:2)],
#'                       w_all = data$treat,
#'                       sl_lib = c("SL.xgboost"),
#'                       dnorm_log = FALSE)
#' }

estimate_gps <- function(cov_mt, w_all, sl_lib, dnorm_log) {

  logger::log_info("Started estimating GPS values ... ")
  t_1 <- proc.time()
  GPS_SL <- SuperLearner::SuperLearner(Y = w_all,
                                       X = as.data.frame(cov_mt),
                                       SL.library = sl_lib)

  GPS_SL_sd <- sd(w_all - GPS_SL$SL.predict)
  GPS_m <- data.frame(GPS = dnorm(w_all,
                                  mean = GPS_SL$SL.predict,
                                  sd = GPS_SL_sd,
                                  log = dnorm_log),
                      e_gps_pred = GPS_SL$SL.predict,
                      e_gps_std = GPS_SL_sd,
                      w = w_all)

  t_2 <- proc.time()
  logger::log_debug("Wall clock time to estimate GPS values:  ",
                    " {t_2[[3]] - t_1[[3]]} s.")

  result <- list()
  class(result) <- "gps"
  result$gps <- GPS_m
  result$used_params <- list(dnorm_log = dnorm_log)

 return(result)
}
