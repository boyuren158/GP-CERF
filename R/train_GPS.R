#' @title
#' Train Model for GPS
#'
#' @description
#' Estimate the conditional mean and sd of exposure level as a function of covariates with
#' xgboost algorithm.
#'
#' @param cov.mt Covariate matrix containing all covariates. Each row is a sample and each
#' column is a covariate.
#' @param w.all A vector of observed exposure levels.
#' @param dnorm_log Logical, if TRUE, probabilities p are given as log(p).
#'
#' @return
#' A data.table that includes:
#'   - a vector of estimated GPS at the observed exposure levels;
#'   - a vector of estimated conditional means of exposure levels when the covariates are fixed
#' at the observed values;
#'   - estimated standard deviation of exposure levels
#' @export
#'
#' @examples
#' mydata <- generate_synthetic_data()
#' GPS_m <- train_GPS(as.matrix(mydata[,c("cf1", "cf2", "cf3", "cf4",
#'                                           "cf5", "cf6")]),
#'                       as.matrix(mydata$treat))
#'
#'
train_GPS <- function(cov.mt, w.all, dnorm_log = FALSE){
  GPS_mod <- xgboost::xgboost(data = cov.mt,
                              label = w.all,
                              nrounds=50,
                              verbose = 0)

  logger::log_info("Started estimating GPS values ... ")
  t_1 <- proc.time()

  e_gps_pred <- predict(GPS_mod,cov.mt)
  e_gps_std <- sd(w.all-e_gps_pred)
  GPS <- c(stats::dnorm(w.all, mean = e_gps_pred, sd = e_gps_std, log = dnorm_log))
  GPS_m <- data.table::data.table(GPS = GPS,
                                  e_gps_pred = e_gps_pred,
                                  e_gps_std = e_gps_std)
  t_2 <- proc.time()
  logger::log_debug("Wall clock time to estimate GPS values:  ",
                    " {t_2[[3]] - t_1[[3]]} s.")

 return(GPS_m)
}
