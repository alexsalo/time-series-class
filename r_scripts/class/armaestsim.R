arma_est <- function(n = c(50, 100, 250, 500), alpha, beta, rvar = 1, 
                   nreps = 500, seed = 0)
#----------------------------------------------------------------------
#  R function to investigate methods of estimation of AR coefficients.
#  NOTE: R stat function "arima" defines the ARMA(p,q) model as
#     x[t] - m = a[1]*x[t-1] + a[2]x[t-2] + ... + a[p]x[t-p] + 
#                + e[t] + b[1]e[t-1] + ... + b[1]e[t-q].
#  We define ARMA(p,q) model as
#     x[t] + a[1]*x[t-1] + a[2]x[t-2] + ... + a[p]x[t-p] = 
#           e[t] + b[1]e[t] + ... + b[q]e[t-q]
#  so it is necessary to take the negative of the estimates of the
#  AR coefficients provided from the arima function.
#
#  INPUT: n = an integer vector of length no more than four containing
#             the lengths of the series.  The default value is
#             n = c(50, 100, 250 500).
#         alpha = a real vector of length no more than three containing 
#                 the AR coefficients.
#         beta = a real vector of length no more than three containing
#                the MA coefficients.
#         rvar = a real scalar containing the error variance.  The
#                default value is rvar = 1.
#         nreps = an integer scalar containing the number of series
#                 to generate (the number of replications).  The 
#                 default value is nreps = 500.
#         seed = an integer scalar containing a seed for the random
#                number generator.
#
#  REQUIRED: The function ar_est requires the following functions that
#            are not a part of base R.
#      Timeslab function armadt(alpha, beta, n, rvar, seed)
#      R stats function ar(x, aic = TRUE, order.max, method)
#----------------------------------------------------------------------
{
  if(length(n) > 4) stop("\nThe length of n can be no more than four.\n")
  if(length(alpha) > 3) stop("\nThe length of alpha can be no more than three.\n")
  if(length(beta) > 3) stop("\nThe length of beta can be no more than three.\n")
#  
  if(length(alpha) == 1 && alpha != 0) {
    zz <- max(Mod(polyroot(c(1,alpha))))
    if(zz <= 1) stop("\nAR characteristic roots not all outside the unit circle.\n")
  }
  if(length(beta) == 1 && beta != 0) {
    zz <- max(Mod(polyroot(c(1,beta))))
    if(zz <= 1) stop("\nMA characteristic roots not all outside the unit circle.\n")
  }
#
#  Initialize objects for simulation:
#
  max.n       <- max(n)
  p           <- length(alpha)
  if(p == 1 && alpha == 0) p <- 0
  q           <- length(beta)
  if(q == 1 && beta == 0) q <- 0
  arma.matrix <- matrix(0, nrow = max.n, ncol = nreps)
  mle.est     <- array(0, dim = c(p+q, length(n), nreps))
  css.est     <- array(0, dim = c(p+q, length(n), nreps))
#  
  for(j in 1:nreps) {
    arma.matrix[1:max.n, j] <- armadt(alpha = alpha, beta = beta, n = max.n, rvar = rvar, seed = seed)$x
    for(i in 1:length(n)) {
      nn <- n[i]
      mle.est[1:(p+q),i,j] <- arima(arma.matrix[1:nn,j], order = c(p,0,q), include.mean = FALSE, method = "ML")$coef
      css.est[1:(p+q),i,j] <- arima(arma.matrix[1:nn,j], order = c(p,0,q), include.mean = FALSE, method = "CSS")$coef
    }
  }
  if(p >= 1 && alpha != 0) {
    mle.est[1:p,,] <- -mle.est[1:p,,]
    css.est[1:p,,] <- -css.est[1:p,,]
  }
#
  return(list(mle.est = mle.est, css.est = css.est))
}