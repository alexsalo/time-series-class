plotpars <- function(alpha = 0, beta = 0, rvar = 1, n = 200, m = 20, Q = 256,
                     ts.main = "Time Series Plot", type = "v",
                     acf.main = "Autocorrelation Function",
                     pacf.main = "Partial Autocorrelation Function",
                     spec.main = "Spectral Density Function")
#----------------------------------------------------------------------
# R function to plot quantities related to a stationary ARMA process.
#
# INPUT: alpha = a real vector of length p containing the 
#                AR coefficients.  The default is alpha = 0.
#        beta = a real vector of length q containing the 
#                MA coefficients. The default is beta = 0.
#        rvar = a real scalar containing the variance of the errors.
#               The default value is rvar = 1.
#        n = an integer scalar indicating the length of the realzation.
#            The default value is n = 200.
#        m = an integer scalar indicating the maximum number of lags
#            for calculating the autocorrelation and partial
#            autocorrelation functions.  The default is m = 20.
#        Q = an integer scalar indicating the number of frequencies on
#            [0, 1] for computing the spectral density function.  The
#            default value is Q = 256.
#        ts.main = a character string containing title for time plot of
#                  realization.  The default value is
#                  ts.main = "Time Series Plot".
#        type = a character indicating what type of plot should be
#               drawn for the acf and pacf.  The default is type = "v". 
#               Possible types are "p" for points, "l" for lines, "b" 
#               for both points and lines, or "v" for vertical lines 
#               with points.
#        acf.main = a character string containing title for the plot of
#                   the autocorrelation function.  The default value is
#                   acf.main = "Autocorrelation Function".
#        pacf.main = a character string containing title for the plot of
#                    the parital autocorrelation function.  The default
#                    value is 
#                    pacf.main = "Partial Autocorrelation Function".
#        spec.main = a character string containing title for the plot of
#                    the spectral density function.  The default value is
#                    spec.main = "Spectral Density Function".
#
#
# RESULT: The function plotpars() creates a 2x2 array of plots.
#          Upper left: a realization of length n of the ARMA process.
#          Upper right: the autocorrelation function.
#          Lower left: the natural logarithn of the standardized 
#                      spectral density function on [0, 0.5].
#          Lower right: the partial autocorrelation function.
#         The function returns a list containing the following objects.
#         p = an integer containing the AR order
#         q = an integer containing the MA order
#         alpha = a real vector of length p containing the  (input)
#                 AR coefficients.
#         beta = a real vector of length q containing the (input) MA
#                coefficients.
#         sigma2 = a real scalar containing the variance of the process.
#         x = a real vector of length n containing a realization of
#             the ARMA(p, q) process.
#         rho = a real vector of length m containing the
#               autocorrelation function.
#         theta = a real vector of length m containing the partial
#                autocorrelation function.
#         f = a real vector of length [Q/2]+1 containing the spectral
#             density function on [0, 0.5].
#
#  Written: 11/18/2015 JLH for exam 2 graphs.
#----------------------------------------------------------------------
{
#  
#  Get ARMA orders p and q.  
#
  p <- length(alpha)
  q <- length(beta)
#
#  Get realization.
#
  x <- armadt(alpha = alpha, beta = beta, n = n, rvar = rvar)
  if(x$error != 0) stop("\n The process is nonstationary.")
  x <- x$x
#
#  Get autocorrelation function.
#
  rho <- armacorr(alpha = alpha, beta = beta, m = m, rvar = rvar)
  if(rho$error != 0) stop("\n The ACF algorithm did not converge.")
  sigma2 <- rho$var
  rho <- rho$corr
#
#  Get partial autocorrelation function.
#
  theta <- armapart(alpha = alpha, beta = beta, rvar = rvar, m = m)
  if(theta$error != 0) stop("\n The PACF algorithm did not converge.")
  theta <- theta$theta
#
#  Get spectral density function.
#
  f <- armasp(alpha = alpha, beta = beta, rvar = rvar, Q = Q)
#
#  Create 2x2 array of plots.
#
  par(mfrow = c(2, 2))
  tsplot(x, main = ts.main)
  plotcorr(r = rho, m = m, corr.type = "acf", type = type, main = acf.main)
  plotcorr(r = theta, m = m, corr.type = "pacf", type = type, main = pacf.main)
  plotsp(f = f, n = Q, div = sigma2, main = spec.main)
#
  return(list(p = p, q = q, sigma2 = sigma2, alpha = alpha, beta = beta,
              x = x, rho = rho, theta = theta, f = f)) 
}