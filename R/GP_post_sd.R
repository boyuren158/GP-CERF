GP_post_sd = function(w, obs.use, param, sigma, e_gps_pred, e_gps_std,
                           kernel.fn = function(x) exp(-x^2)){
  # param[1]: alpha, param[2]: beta, param[3]: gamma
  n = nrow(obs.use)
  GPS.new = stats::dnorm(w, mean = e_gps_pred, sd = e_gps_std, log = T)
  
  obs.new = cbind( w/sqrt(param[1]), GPS.new/sqrt(param[2]) )
  obs.all = rbind(obs.new, obs.use)
  Sigma.all = (param[3]*kernel.fn(stats::dist(obs.all)) + diag(n*2))*sigma^2
  Sigma.within.w = Sigma.all[1:n, 1:n]
  Sigma.cross = Sigma.all[1:n, -(1:n)]
  Sigma.within.obs = Sigma.all[-(1:n), -(1:n)]
  Sigma.cond = Sigma.within.w - Sigma.cross%*%chol2inv(Sigma.within.obs)%*%t(Sigma.cross)
  sqrt(rep(1/n,n)%*%Sigma.cond%*%rep(1/n,n) + sigma^2)
}