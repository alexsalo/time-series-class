ar_est <- function(n = c(50, 100, 250, 500), alpha, rvar = 1, 
                   nreps = 500, seed = 0)
#----------------------------------------------------------------------
#  R function to investigate methods of estimation of AR coefficients.
#  NOTE: R stat function "ar" defines the autoregressive model as
#     x[t] - m = a[1]*x[t-1] + a[2]x[t-2] + ... + a[p]x[t-p] + e[t].
#  We define AR model as
#     x[t] + a[1]*x[t-1] + a[2]x[t-2] + ... + a[p]x[t-p] = e[t],
#  so it is necessary to take the negative of the estimates provided
#  from the ar function.
#
#  INPUT: n = an integer vector of length no more than four containing
#             the lengths of the series.  The default value is
#             n = c(50, 100, 250 500).
#         alpha = a real vector of length no more than three containing 
#                 the AR coefficients.
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
#      Timeslab function ardt(alpha, n, rvar, seed)
#      R stats function ar(x, aic = TRUE, order.max, method)
#----------------------------------------------------------------------
{
  if(length(n) > 4) stop("\nThe length of n can be no more than four.\n")
  if(length(alpha) > 3) stop("\nThe length of alpha can be no more than three.\n")
#  
  zz <- max(Mod(polyroot(c(1,alpha))))
  if(zz <= 1) stop("\nThe AR model must be stationary.\n")
#
#  Initialize objects for simulation:
#
  max.n     <- max(n)
  p         <- length(alpha)
  ar.matrix <- matrix(0, nrow = max.n, ncol = nreps)
  yw.est    <- array(0, dim = c(p, length(n), nreps))
  ols.est   <- array(0, dim = c(p, length(n), nreps))
  mle.est   <- array(0, dim = c(p, length(n), nreps))
  burg.est  <- array(0, dim = c(p, length(n), nreps))
#  
  for(j in 1:nreps) {
    ar.matrix[1:max.n, j] <- ardt(alpha = alpha, n = max.n, rvar = rvar, seed = seed)$x
    for(i in 1:length(n)) {
      nn <- n[i]
      yw.est[1:p,i,j]   <- -ar(ar.matrix[1:nn,j], AIC = FALSE, order.max = p, method = "yule-walker")$ar
      ols.est[1:p,i,j]  <- -ar(ar.matrix[1:nn,j], AIC = FALSE, order.max = p, method = "ols")$ar
      mle.est[1:p,i,j]  <- -ar(ar.matrix[1:nn,j], AIC = FALSE, order.max = p, method = "mle")$ar
      burg.est[1:p,i,j] <- -ar(ar.matrix[1:nn,j], AIC = FALSE, order.max = p, method = "burg")$ar
    }
  }
#
  return(list(yw.est = yw.est, ols.est = ols.est, mle.est = mle.est, 
              burg.est = burg.est))
}