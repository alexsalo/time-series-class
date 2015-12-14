arcorr <- function(alpha, m = length(alpha) + 10, rvar = 1)
#-----------------------------------------------------------------------
#  R function to evaluate autocorrelation function and variance for an  
#  AR process with coefficients alpha and error variance rvar.
#
#  INPUT: alpha = a real vector containing AR coefficients.
#         m = an integer scalar containing maximum autocorrelation lag.
#             The default value is length(alpha) + 10.
#         rvar = a real scalar containing error variance.  The
#                default value is rvar = 1.
#
#  VALUE: The function arcorr returns a list containing the following
#         objects:
#         var = a real scalar containing the variance of the AR 
#               process.
#         corr = a real vector of length m containing the
#                autocorrelation function of the AR process.
#         error = an integer scalar that is an error indicator.  If 
#                 error = 0, there was no error.  If error = 1, then
#                 the value of p or q was negative.  If error = 2,
#                 then the value of m was less than max(p,q)+1.
#                 If error = 3, then the AR model was not stationary.
#
#  FORTRAN required : SUBROUTINE wilson (within armacorr function)
#-----------------------------------------------------------------------
{
  return(armacorr(alpha = alpha, beta = 0, m, rvar = rvar))
}
#
#
#
ardt <- function(alpha, n, rvar = 1, seed = 0)
#-----------------------------------------------------------------------
#  R function to simulate a realization of length n from an AR process 
#  with coefficients in alpha, error variance rvar.
#
#  INPUT: alpha = a real vector containing the coefficients of the AR
#                 process.
#         n = an integer scalar indicating the length of the desired
#             realization.
#         rvar = a real scalar containing the variance of the errors.
#                The default value is rvar = 1.
#         seed = a nonnegative integer for the random number generator.
#                If seed = 0 (the default) then the output seed from
#                the last call is used.  "seed" is an optional argument.
#
#  VALUE: The function ardt returns a list containing the following
#         two objects.
#         x = a real vector of length n containing the desired 
#             realization.
#         error = an integer scalar that is an error indicator.  If 
#                 error = 0, there was no error.  If error = 1, then
#                 the value of p or q was negative.  If error = 2,
#                 then the value of m was less than max(p,q)+1.
#                 If error = 3, then the AR model was not stationary.
#
#  FORTRAN required: SUBROUTINE wilson (within arcorr function)
#-----------------------------------------------------------------------
{
  x  <- rep(0, length = n)
  p  <- length(alpha)
  z  <- arcorr(alpha, p, rvar)
 
  if(z$error == 0) {
    e <- c(corrdt(z$corr, z$var, p, seed), 
            sqrt(rvar)*wn(n-p, dist = "normal", seed = seed))
    x <- diffeq(alpha, n, e, p)
  }
  return(list(x = x, error = z$error))
}
#
#
#
ar.eigen <- function(alpha, n = length(alpha)+1, rvar = 1)
#-----------------------------------------------------------------------
#  An R function to compute the n eigenvalues of covariance matrix and
#  sorted spectral density function of an AR(p) process with 
#  coefficients as specified by alpha.
#
#  INPUT: alpha = a real vector containing the coefficients of the AR
#                 process
#         n = an integer scalar containing the number of eigenvalues
#             and spectral density values to compute.  The value
#             of n must be at least length(alpha) + 1.  The default
#             value is n = length(alpha) + 1.  (Note that the maximum
#             lag will be n - 1, and the dimension of the Toeplitz
#             matrix will by n x n.)
#         rvar = a real scalar containing the variance of the errors.
#                The default value is rvar = 1.
#
#  VALUES: The function ar.eigen returns a list containing the following
#          objects.
#          eigenvals = a vector of length n containing the eigenvalues
#                      of the AR(p) process if the process is 
#                      covariance stationary.  If not, eigenvals is NULL.
#          sorted.f = a vector of length n containing the sorted values
#                     of the spectral density function if the process
#                     is second-order stationary.  If the process is
#                     not stationary, sorted.f = NULL.
#          error = an integer indicator.  If error = 0, then the
#                  process is second-order stationary.  Otherwise
#                  error = 1.
#-----------------------------------------------------------------------
{
  if(n < length(alpha)) n <- length(alpha) + 1
  m <- n - 1

  cc <- arcorr(alpha, m, rvar)
  error <- cc$error
  if(error == 0) {
    R0 <- cc$var
    R  <- cc$corr * R0
    evals <- eigen(toepl(R, R0), only.values = TRUE, symmetric = TRUE)$values

    f1 <- arsp(alpha, rvar, n)
    f  <- sort(c(f1[1:(length(f1)-1)], f1[1:(length(f1)-1)]), decreasing = TRUE)

  }
  if(error != 0) {
    evals <- NULL
    f     <- NULL
  }

  return(list(eigenvals = evals, sorted.f = f, error = cc$error))

}
#
#
#
armacorr <- function(alpha, beta, m = (length(alpha) + length(beta) + 10), 
                     rvar = 1)
#-----------------------------------------------------------------------
#  R function to evaluate autocorrelation function and variance for an  
#  ARMA process with AR coefficients alpha, MA coefficients beta, and 
#  error variance rvar.
#
#  INPUT: alpha = a real vector containing the AR coefficients.
#         beta = a real vector containing the MA coefficients.
#         m = an integer scalar containing maximum autocorrelation lag.
#             The default value is m = length(alpha) + length(beta) + 10.
#         rvar = a real scalar containing error variance.
#
#  VALUE: The function arcorr returns a list containing the following
#         objects:
#         var = a real scalar containing the variance of the ARMA 
#               process.
#         corr = a real vector of length m containing the
#                autocorrelation function of the ARMA process.
#         error = an integer scalar that is an error indicator.  If 
#                 error = 0, there was no error.  If error = 1, then
#                 the value of p or q was negative.  If error = 2,
#                 then the value of m was less than max(p,q)+1.
#                 If error = 3, then the ARMA model was not stationary.
#
#  FORTRAN required: SUBROUTINE wilson
#-----------------------------------------------------------------------
{
  p <- length(alpha)
  q <- length(beta)
  if(p == 0) alpha <- 0
  if(q == 0) beta  <- 0

  ier   <- 0

  z <- .Fortran("wilson",
                as.double(-alpha),
                as.integer(p),
                as.double(-beta),
                as.integer(q),
                as.integer(m+1),
          acf = as.double(rep(0, length = (m+1))),
          ier = as.integer(ier))
#
  error <- z$ier
  tsvar <- rvar*z$acf[1]
  if(m == 0 || error != 0) 
    corr <- NULL
  else 
    corr <- z$acf[2:(m+1)] / tsvar

  return(list(var = tsvar, corr = corr, error = error))
}
#
#
#
armadt <- function(alpha, beta, n, rvar = 1, seed = 0)
#-----------------------------------------------------------------------
#  R function to simulate a realization of length n from an ARMA process
#  having AR coefficients alpha, MA coefficients beta, and error
#  variance rvar.  
#
#  INPUT: alpha = a real vector containing the coefficients of the AR
#                 part of the model.
#         beta = a real vector containing the coefficients of the MA 
#                part of the model.
#         n = an integer scalar specifying the length of the desired
#             realization.
#         rvar = a real scalar containing the variance of the error
#                terms.  The default value is rvar = 1.
#         seed = a nonnegative integer for the random number generator.
#                If seed = 0 (the default) then the output seed from
#                the last call is used.  "seed" is an optional argument.
#
#  VALUE: The function armadt returns a list containing the following
#         objects.
#         x = a real vector of length n containing the desired ARMA
#             realization.  If the ARMA process is nonstationary, then
#             x will be empty.
#         error = an integer indicating whether the ARMA process is 
#                 stationary.  error = 0 indicates the process is
#                 stationary.  Any other value indicates it is not.
#
#  FORTRAN required: diffeq
#-----------------------------------------------------------------------
{
  error <- 0
  p <- length(alpha)
  q <- length(beta)

  if(p == 0 && q == 0) x <- rnorm(n, sqrt(rvar))
  if(p == 0) x <- madt(beta, rvar, n, seed)
  if(q == 0) {
     xx <- ardt(alpha, n, rvar, seed)
     x  <- xx$x
     error <- xx$error
  }
  if(p != 0 && q != 0) {
    xx <- armacorr(alpha = alpha, beta = beta, m = (p+q), rvar = rvar)
    if(xx$error != 0) {
      x <- c()
      error <- xx$error
    }
    else { 
      e <- c(corrdt(xx$corr, xx$var, n = p, seed = seed), 
                madt(beta = beta, n = (n - p), rvar = rvar,
                seed = seed))
      x <- diffeq(alpha = alpha, n = n, e = e, p = p)
    }
  }

  return(list(x = x, error = error))
}
#
#
#
armapart <- function(alpha, beta, rvar = 1, m = max(length(alpha), length(beta))+10)
#-----------------------------------------------------------------------
#  R function to compute ARMA partial partial autocorreltions.
#
#  INPUT: alpha = a real vector containing the coefficients for the
#                 AR part of the model.
#         beta = a real vector containing the coefficients for the
#                MA part of the model.
#         rvar = a real scalar containing the error variance.  The
#                default value is rvar = 1.
#         m = an integer scalar specifying the number of partial
#             autocorrelations to return.  The default value is 
#             m = max(p, q) + 10.
#
#  VALUE: The function armapart returns a list containing the following
#         objects:
#         theta = a real vector of length m containing the partial
#                 autocorrelations of the ARMA process.  If error is
#                 not zero, then theta is an empty vector.
#         error = an integer scalar indicating if there was no error
#                 (error = 0) or if the specified model was 
#                 nonstationary (error not 0).
#------------------------------------------------------------------------
{
   theta <- c()
   z <- armacorr(alpha, beta, m, rvar)
   error <- z$error
   if(error == 0) {
    z <- corrar(z$corr, z$var)
    theta <- arpart(z$alpha)$theta
   }
   return(list(theta = theta, error=error))
}
#
#
#
arpart <- function(alpha)
#-----------------------------------------------------------------------
#  R function to find partial autocorrelation function of an AR process
#  having coefficients alpha.
#
#  INPUT: alpha = a real vector containing AR coefficients.
#
#  VALUE: The function arpart returns list containing the following 
#         objects.
#         theta = a real vector of length p containing the partial 
#                 autocorrelation function for lags 1, ..., p.
#         error = an integer scalar indicating if the model was
#                 stationary (error = 0) or not (error not 0).
#-----------------------------------------------------------------------
{
  p <- length(alpha)
  theta <- rep(0, length = p)
  if(p != 0) {
    theta  <- lambda <- -alpha[1:p]
    for(j in p:1)
    {
      if(abs(theta[j]) >= 1) return(list(theta = theta, ier = j))
      if(j == 1) return(list(theta = theta, ier = 0))
      lambda <- lambda[1:(j-1)]
      lambda <- (lambda + theta[j] * rev(lambda)) / (1 - theta[j]^2)
      theta[(j-1)] <- lambda[(j-1)]
   }

  return(list(theta=theta,ier=0))
}
}
#
#
#
arsp <- function(alpha, rvar = 1, Q = 256)
#-----------------------------------------------------------------------
#  R function to compute the spectral density function at Q frequencies
#  in [0, 1] for an AR(p) process with coefficients alpha and error 
#  rvar.
#
#  INPUT: alpha = a real vector containing the coefficients of the AR(p)
#                 process.
#         rvar = a real scalar containing the variance of the errors.
#                The default value is rvar = 1.
#         Q = an integer scalar indicating the number of frequencies
#             in the interval [0,1] at which to compute the spectral 
#             density function.  The default value is Q = 256.
#
#  VALUE: The function arsp returns a real vector of length [Q/2] + 1
#         containing the spectral density function of the AR(p) 
#         process.
#-----------------------------------------------------------------------
{
  return(armasp(alpha = alpha, beta = NULL, rvar = rvar, Q = Q))
}
#
#
#
armasp <- function(alpha, beta, rvar = 1, Q = 256)
#-----------------------------------------------------------------------
#  R function to compute the spectral density function at Q frequencies
#  in [0, 1] for an ARMA(p, q) process with AR coefficients alpha,
#  MA coefficients beta, and error variance rvar.
#
#  INPUT: alpha = a real vector containing the coefficients of the AR
#                 part of the ARMA process.
#         beta = a real vector containing the coefficients of the MA
#                part of the ARMA process.
#         rvar = a real scalar containing the variance of the errors.
#                The default value is rvar = 1.
#         Q = an integer scalar indicating the number of frequencies
#             in the interval [0, 1] at which to compute the spectral 
#             density function.  The default value is Q = 256.
#
#  VALUE: The function armasp returns a real vector of length [Q/2]+1
#         containing the spectral density function of the ARMA(p, q) 
#         process.
#-----------------------------------------------------------------------
{
  f1 <- f1 <- rep(1, length = Q)

  p <- length(alpha)
  q <- length(beta)

  if(q != 0) f1 <- Re(Mod(fft(c(1, beta, rep(0, length = (Q - q - 1)))))^2)
  if(p != 0) f2 <- Re(Mod(fft(c(1, alpha, rep(0, length = (Q - p - 1)))))^2)

  f <- rvar * f1 / f2
  f <- f[1:(floor(Q/2)+1)]

  return(f)
}

#
#
#  
barttest <- function(F, n) 
#-----------------------------------------------------------------------
#  R function to compute values of the Bartlett test for white noise
#  based on the cumulative periodogram of a data set of length n.
#
#  INPUT: F = a real vector containing the cumulative periodogram of
#             the time series or a real scalar containing the Bartlett
#             test statistic for which a p-value is desired.
#         n = an integer scalar containing the length of the original
#             series.
#
#  VALUE: The function barttest returns a list containing the following
#         objects.
#         B = a real scalar containing the Bartlett test statistic.
#         pval = a real scalar containing the p-value of the Bartlett
#                statistic for testing white noise.
#-----------------------------------------------------------------------
{
  n1 <- length(F)
  if(n1 == 1) B <- F else
    B  <- sqrt(n/2)*max(abs(F[2:n1] - ((2:n1)/n1)))

  pval <- -2*sum((-1)^(1:50)*exp(-2*B^2*(1:50)^2))

  return(list(B = B, pval = pval))
}
#
#
#
clip <- function(x, low = -1.e20, up = 1.e20)
#-----------------------------------------------------------------------
#  R function to clip an array x so that all values are between 
#  low and up.
#
#  INPUT: x = a real vector (or scalar) containing the data to be 
#             clipped.
#         low = a real scalar containing the lower bound.  The default
#               is -1 x 10^{20}.
#         up = a real scalar containing the upper bound.  The default
#              is 1 x 10^{20}.
#
#  VALUE: The function clip returns an real vector of the same length as
#         x with values less than low replaced by low and values greater
#         than up replaced by up.
#------------------------------------------------------------------------
{
   ind <- c(1:length(x))
   x   <- replace(x, ind[x < low], low)
   x   <- replace(x, ind[x > up], up)
   return(x)
}
#
#
#
corr <- function(x, m = 0)
#-----------------------------------------------------------------------
#  R function to calculate the sample autocorrelation function (acf)
#  lags 1 through m via 2 FFTs.
#
#  INPUT: x = a real vector containing the time series.
#         m = an integer scalar containing the number of lags
#             at which to find the sample acf.  The default is
#             m = 0, in which case only the sample time series
#             variance is returned.
#
#  VALUE: If m = 0, the corr function returns a real scalar containing
#         the time series sample variance.  If m is not equal to 0,
#         the function returns a list containing the following objects:
#         corr = a real vector of length m containing the sample acf,
#         var = a real scalar containing the time series variance.
#-----------------------------------------------------------------------
{
  z <- x - mean(x)
  if(m == 0) xx <- tsvar(z) 
  else {
    n   <- length(x)
    Q   <- 2^(trunc(log(n + m, base = 2) + 0.05) + 1)
    if(Q < (n + m)) Q <- 2*Q
    z   <- c(z, rep(0, length = (Q-n)))
    fz  <- Re( fft( Mod(fft(z))^2 / n, inverse = TRUE)) / Q
#
    xx  <- list(corr = fz[2:(m + 1)]/fz[1], var = fz[1])
  }
  return(xx)
}
#
#
#
corr1 <- function(x, m = 0)
#-----------------------------------------------------------------------
#  R function to calculate the sample autocorrelation function (acf)
#  lags 1 through m via the convolution method.
#
#  INPUT: x = a real vector containing the time series.
#         m = an integer scalar containing the number of lags
#             at which to find the sample acf.  The default is
#             m = 0, in which case only the sample time series
#             variance is returned.
#
#  VALUE: If m = 0, the corr function returns a real scalar containing
#         the time series sample variance.  If m is not equal to 0,
#         the function returns a list containing the following objects:
#         corr = a real vector of length m containing the sample acf,
#         var = a real scalar containing the time series variance.
#-----------------------------------------------------------------------
{
  var <- tsvar(x)
  if(m == 0) xx <- var 
  else {
    n   <- length(x)
    x   <- c( x - mean(x), rep(0, length = m) )
    y   <- x
    rho <- rep(0, length = m)
    for(i in 1:m) 
      rho[i] <- dot( x, (y <- crlag(y))) / (n * var)
#
    xx <- list(corr = rho, var = var)
  }
  return(xx)
}
#
#
#
corr2 <- function(x, y, m = 0)
#-----------------------------------------------------------------------
#  Function to compute sample cross-correlation between two series x 
#  and y.
#
#  INPUT: x, y = two real vectors containing time series of equal 
#                lengths between which the sample cross-correlation 
#                is desired.
#         m = an integer scalar containing the number of lags at which
#             to compute the sample cross-correlation function. The
#             default is m = 0.
#
#  VALUE: If m = 0, then the function corr2 returns a list containing
#         three objects:
#         var.x = a real scalar containing the sample times series
#                 variance of x
#         var.y = a real scalar containing the sample times series
#                 variance of y
#         rho = a real scalar containing the sample correlation
#               coefficient between x and y.
#         If m is not equal to zero, then the function corr2 returns 
#         a list containing the following four objects.
#         var.x = a real scalar containing the sample times series
#                 variance of x
#         var.y = a real scalar containing the sample times series
#                 variance of y
#         lags = an integer vector containing the lags at which
#                the sample cross-correlation was computed (-m:m)
#         cross.corr = a real scalar vector of length 2m + 1 
#                containing the sample cross-correlation between
#                the two series x and y.
#-----------------------------------------------------------------------
{
  if(length(x) != length(y)) 
    stop("\nThe lengths of the two series must be equal.\n")

  rho  <- cor(x, y)
  varx <- tsvar(x)
  vary <- tsvar(y)

  if(m == 0) xx <- list(var.x = varx, var.y = vary, corr = rho) else {
    n <- length(x)
    x <- x1 <- c(x - mean(x), rep(0, length = m))
    y <- c(y - mean(y), rep(0, length = m))
    rho.xy <- rep(0, length = (2*m + 1))

    div <- length(x)*sqrt(varx*vary)
    rho.xy[m+1] <- rho

    for(i in 1:m) rho.xy[m+1+i] <- dot(y, x1 <- crlag(x1)) / div
    for(i in m:1) rho.xy[i] <- dot(x, y <- crlag(y)) / div

    lags <- -m:m
#
    xx <- list(var.x = varx, var.y = vary, lags = lags, cross.corr = rho.xy)
  }
  return(xx)
}
#
#
#
corrar  <- function(rho, r0)
#-----------------------------------------------------------------------
#  R function to find AR parameters given the autocorrelation function.
#
#  INPUT: rho = a real vector of length p containing the autocorrelation
#               autocorrelation function.
#         r0 = a real scalar containing the variance of the series.
#
#  VALUE: The function corrar returns a list containing the following
#         objects:
#         alpha = a real vector of length p containing the AR
#                 coefficients.
#         rvar = a real scalar containing the variance of the error
#                series.
#----------------------------------------------------------------------
{
  p     <- length(rho)
  alpha <- rep(0, length = p)
  rvar  <- 0

  z <- .Fortran("corrar", 
                 as.double(rho), 
                 as.double(r0),
                 as.integer(p),
         alpha = as.double(alpha),
          rvar = as.double(rvar))

   return(list(alpha = z$alpha, rvar = z$rvar))
}
#
#
#
corrdt <- function(rho, r0, n, seed = 0)
#-----------------------------------------------------------------------
#  R function to simulation a Gaussian time series realization of 
#  length n having an autocorrelation function as specified in rho and 
#  variance r0.
#
#  INPUT: rho = a real vector containing the values of the
#               autocorrelation function.
#         r0 = a real scalar containing the variance of the series.
#         n = an integer scalar containing the length of the series.
#         seed = a nonnegative integer for the random number generator.
#                If seed = 0 (the default) then the output seed from
#                the last call is used.  "seed" is an optional argument.
#
#  VALUE: The function corrdt returns a real vector of length n 
#         containing the realization of the process.
#----------------------------------------------------------------------
{
  for(v in 1:length(rho)) {
    if(abs(rho[v]) > 1) 
      stop("\n All values of rho must be between -1 and 1.\n")
  }
  x <- rep(0, length = n)
  if(n != 0) 
    x <- t(chol(toepl(rho*r0, r0, n))) %*% wn(n, dist = "normal", seed = seed)

  return(x)
}
#
#
#
crlag <- function(x)
#-----------------------------------------------------------------------
#  R function to apply the circular shift operator to an array x. 
#
#  INPUT: x = a real vector.
#
#  VALUE: The function crlag returns a real vector that is the
#         result of applying the circular shift operator to x.
#-----------------------------------------------------------------------
{
   return(c(x[length(x)], x[1:(length(x) - 1)]))
}
#
#
#
cum <- function(x, iopt = 1)
#-----------------------------------------------------------------------
#  R function to return the cumulative sum of x if iopt = 1 or
#  cumulative averages is iopt = 2.
#
#  INPUT: x = a real vector containing the data.
#         iopt = an integer scalar indicating whether to compute
#                the cumulative sum (iopt = 1, the default), or
#                cumulative average (iopt = 2).  If iopt is not
#                1 or 2, then iopt = 1 is assumed.
#
#  VALUE: The function cum returns a real scalar containing the
#         cumulative sum if iopt = 1 or the cumulative average 
#         if iopt = 2.
#-----------------------------------------------------------------------
{
  if(iopt != 1 && iopt != 2) iopt <- 1

  if(iopt == 1) xx <- cumsum(x)
  if(iopt == 2) xx <- cumsum(x)/(1:length(x))
#
  return(xx)
}
#
#
#
cumsp <- function(f)
#-----------------------------------------------------------------------
#  R function to compute the cumulative spectral density function.
#
#  INPUT: f = a real vector containing the spectral density function.
#
#  VALUE: The function cumsp returns a vector containing the
#         cumulative spectral density function.
#-----------------------------------------------------------------------
{
  return(cum(f) / sum(f))
}
#
#
#
descplot <- function(x, m, alpha = 0.05)
#-----------------------------------------------------------------------
#  R function to create plots of univariate descriptive statistics.
#
#  INPUT: x = a real vector containing the time series
#         m = an integer scalar indicating the number of lags for
#             computing the sample autocorrelation function, sample
#             partial autocorrelation functions, and standardized
#             residual variances.
#         alpha = a real scalar between 0 and 0.2 indicating the level
#                 of significance for white noise tests.
#  VALUE: The function descplot creates plots of univariate descriptive
#         statistics of the time series in x.
#-----------------------------------------------------------------------
{
  if(m < 5) stop("\n The value of m must be greater than 5.\n")
  if(alpha < 0 || alpha > 0.20)
    stop("\nThe value of alpha must be between 0 and 0.2.\n")

  z <- describe(x, m)
  n <- length(x)

  par(mfrow = c(3, 2), oma = rep(0, length = 4), mar = c(5, 4, 2, 2) + 0.1)
  plot(x, type = "l", xlab = "t", ylab = expression(x[t]), 
       main = "Time Plot")
  plot(z$corr, ylim = c(-1,1), type = "l",
       xlab = "v", ylab = expression(hat(rho[v])), 
       main = "Correlogram")
  lim <- qnorm((1+(1-alpha)^(1/m))/2)/sqrt(n)
  lines(c(0, m), c(lim, lim), lty = 2)
  lines(c(0, m), c(-lim, -lim), lty = 2)
  lines(c(0, m), c(0, 0))
#
  plot(z$part, ylim = c(-1, 1), type = "l",
       main = "Partial Correlogram",
       xlab = "v", ylab = expression(hat(theta[v])))
  lines(c(0, m), c(0, 0))
#
  plot(z$srvar, ylim = c(0, 1), pch = 20,
       main ="Standardized Residual Variances",
       xlab = "v", ylab = expression(hat(sigma^2)[v]))
  plot(freqs(n), log(stdf(z$per, tsvar(x), exp(-6), exp(6))),
       type = "l", xlab = expression(omega), ylab = "", ylim = c(-6, 6),
       main="Natural Log of Standardized Periodogram")
  plot(freqs(n), z$cper, ylim = c(0, 1), type = "l",
       xlab = expression(omega), ylab = expression(hat(F)*(omega)),
       main = "Cumulative Periodogram")
#
  fb <- rep(0, length = 501)
  b  <- seq(0.4, 2, length = 501)
  for(i in 1:501) fb[i] <- 1 - barttest(b[i], n)$pval
  inds <- which(round(fb, 3) == round(1 - alpha, 3))
  if(length(inds) == 1) cc <- b[inds] else cc <- median(b[inds])
  xx <- seq(0, 0.5, length = n)
  yy.upper <- 2*xx + cc/sqrt(length(z$cper))
  yy.lower <- 2*xx - cc/sqrt(length(z$cper))

  lines(c(0, 0.5), c(0, 1))
  lines(xx, yy.upper, lty = 2)
  lines(xx, yy.lower, lty = 2)
#
  return(invisible())
}
#
#
#
describe <- function(x, m)
#-----------------------------------------------------------------------
#  R function to compute descriptive statistics for a time series x.
#
#  INPUT: x = a real vector containing the time series.
#         m = an integer scalar containing the number of lags for
#             computing the sample autocorrelation, partial
#             autocorrelation, and standardized residual variances.
#
#  VALUE: The function describe returns a list containing the following
#         objects.
#         mean = the sample mean of the time series in c.
#         var = the time series sample variance of x.
#         per = the periodogram of x.
#         cper = the cumulative periodogram of x.
#         corr = the sample autocorrelation function at lags 1, ..., m.
#         part = the sample partial autocorrelation function for 
#                lags 1, ..., m.
#         srvar = the standardized residual variances at lags 1, ..., m.
#-----------------------------------------------------------------------
{
   per <-  perdgm(x)
   part <- parcorr(x,m)
#
   list(mean = mean(x), var = tsvar(x), per = per, cper = cumsp(per), 
        corr = corr(x,m)$corr, part = part, srvar = srvars(x, m))
}
#
#  
#
diffeq <- function(alpha, n, e, p = length(alpha)) 
#-----------------------------------------------------------------------
#  R function to find realization from a difference equation with 
#  forcing term e.
#
#  INPUT: alpha = real vector containing coefficients of the 
#                 difference equation.
#         n = Integer > length(alpha) containing the length of the
#             realization from the difference equation.
#         e = Array of length n containing the values of the forcing
#             term.
#         p = an integer scalar.  The default value is 
#             p = length(alpha). If p == 0, then alpha is set to zero.
#
#  VALUE: The function diffeq returns an vector of length n containing
#         the n values of the series based on the difference
#         equation.
#
#  FORTRAN required: SUBROUTINE diffeq
#-----------------------------------------------------------------------
{
  if(p < 0 || (as.integer(p) - p) != 0) 
    stop("\n The value of p must be a non-negative integer.\n")
  if(p == 0) alpha <- c(0)
  if(p == 1 && alpha == 0) p <- 0

  z <- .Fortran("diffeq", 
                 as.double(alpha), 
                 as.integer(p), 
                 as.integer(n),
                 as.double(e), 
             x = as.double(rep(0, length = n)))
#
  return(z$x)
}
#
#
#
dirichlet <- function(m, n = 501){
#------------------------------------------------------------------------
#  R function to compute the Dirichlet kernel at value of m.
#
#  INPUT: m = a positive integer.
#         n = a positive integer indicating at how many points the
#             kernel should be computed.  The default is n = 501.
#
#  VALUE: The function dirichlet returns a vector of length n containing
#         values of the Dirichlet kernel for n frequencies from
#         -0.5 to 0.5.
#------------------------------------------------------------------------
  if(n <= 0 || (floor(n)-n) != 0) 
    stop("\n The value of n must be an integer greater than one.\n")
  if(min(m) <= 0 || min(floor(m) - m) != 0)
    stop("\n All values of m must be positive integers.\n")

  omega <- seq(-0.5, 0.5, length = n)

  return(sin(2*pi*omega*(m+0.5))/sin(pi*omega))
}
#
#
#
divpoly <- function(num, den, n)
#-----------------------------------------------------------------------
#  R function to find first n coefficients of the quotient of two 
#  polynomials g(z)/h(z).
#
#  INPUT: num = a real vector containing the coefficients of the 
#               polynomial in the numerator.
#         den = a real vector containing the coefficients of the
#               polynomial in the denominator.
#         n = an integer scalar inidicating the number of coefficients
#             of the quotient polynomial to be returned.
#
#  VALUE: The function divpoly returns a vector containing the first 
#         n coefficients of the quotient polynomial g(z)/h(z).
#
#  FORTRAN required: SUBROUTINE mxma
#-----------------------------------------------------------------------
{
  q <- length(num)
  p <- length(den)
  ratio <- rep(0, length(n))

  if((p+q) != 0) {
    if(p == 0) alpha <- 0 else alpha <- den
    if(q == 0) beta <- 0 else beta <- num
    z <- .Fortran("mxma",
                  as.double(alpha),
                  as.double(beta),
                  as.integer(p),
                  as.integer(q),
                  as.integer(n),
          ratio = as.double(ratio <- rep(0, length=n)))

    ratio <- z$ratio
  }
  return(ratio)
}
#
#
#
divsds <- function(x, d = 1, option = 1, sds = rep(1, length.out = d))
#----------------------------------------------------------------------
#  R function to divide a time series data x by the seasonal standard
#  deviations.
#
#  INPUT: x = a real vector containing the time series.  The length
#             of x must be a multiple of d.
#         d = an integer scalar indicating the length of the cycle;
#             for example, for monthly data across a year, d = 12.
#             The default value of d = 1.  If d = 1, the function
#             computes the sample standard deviation of the series.
#         option = an integer scalar indicating what operation should
#             be performed.  The default value of option is option = 1.
#         sds = a vector of length d containing the seasonal standard
#               deviations. The default value of sds is 
#               rep(0, times = d). sds is ignored if option = 0 or 1.
#
#  The result of divsds is dependent upon the value of the option
#  argument.
#  If option = 0: the function will return a list containing 
#       x = a real vector containing the original series, and
#       sds = a vector of length d containing the seasonal standard
#             deviations.
#     option = 1: the function will return a list containing
#       x = a real vector containing the original series with the
#             divided by the seasonal standard deviations.
#       xbar = a vector of length d containing the seasonal standard
#              deviations.
#     option = 2: the function will return the series x multiplied by
#       the seasonal seasonal standard deviations provided in sds.
#----------------------------------------------------------------------
{
  n <- length(x)
  y <- x
  if(d == 1) {
    sds <- sqrt(tsvar(x))    
    if(option == 0) return(sds)
    if(option == 1) return(list(x = x/sds, sds = sds))
    if(option == 2) return(sds[1]*x)
  }
#  
  if(option == 2) {
    for(i in 1:d) {
      inds <- seq(from = i, to = n, by = d)
      x[inds] <- sds[i]*x[inds]}
    return(x)
  }
#  
  for(i in 1:d) {
    inds <- seq(from = i,to = n, by = d)
    sds[i] <- sqrt(tsvar(x[inds]))
    if(option == 1) x[inds] <- x[inds]/sds[i]
  }
  return(list(x = x, sds = sds))
}
#
#
#
dot <- function(x, y)
#-----------------------------------------------------------------------
#  R function to find the inner product of two vectors, x and y.
#
#  INPUT: x, y = two real vectors of equal length for which the dot
#                product is desired.
#
#  VALUE: The function dot returns a real scalar equal to the dot
#          product of the real vectors x and y.
#-----------------------------------------------------------------------
{
  return(sum(x * y))
}
#
#
#
expsm <- function(x, first = 0, last = 1, k = 101)
#-----------------------------------------------------------------------
#  R function to find value of alpha in exponential smoothing, and then
#  use the exponential smoother to find fitted values, and superimpose
#  those as the same plot as the original data.
#
#  INPUT: x = a real vector containing the data to be smoothed.
#         first = a real scalar containing the minimum value of 
#                 alpha to be considered.  The default is first = 0.
#         last = a real scalar containing the maximum value of alpha
#                to be considered.  The default is last = 1.
#         k = an integer scalar containing the number of values of alpha
#             to consider.  The default is 101.
#
#  VALUE: The function expsm creates a set of two plots.  The first
#         is a plot of values of the sum of squared errors for varying
#         values of the smoothing constant (a) versus a.  The second
#         plot is a plot of the original series (line) with the 
#         fitted values using exponential smoothing with value alpha
#         = min_a (SS(a)) superimposed as points.  The function also
#         returns a list containing the following objects.
#         alpha = a real scalar between first and last containing the
#                 value of the smoothing constant the minimizes the 
#                 sum of squared errors.
#         min.SS = a real scalar containing the minimum value of the
#                  sum of squared errors for all considered values
#                  of smoothing constant.
#         xhat = a real vector of length n containing the fitted
#                values using exponential smoothing, based on the value
#                alpha.
#         avalues = a vector of length k containing all considered
#                   values of smoothing constant.
#         all.SS = a vector of length k containing the error sum of
#                  squares for all a values.
#
#  FORTRAN subroutine: diffeq (through diffeq function)
#
#  Written for R: 10/29/2010 JLH
#-----------------------------------------------------------------------
{
  if(last < first) {
    cc <- last
    last <- first
    first <- cc
  }
  if(first < 0 || last > 1) 
    stop("\n The value of alpha must be between 0 and 1.\n")

  n     <- length(x)
  SS    <- rep(0, length = k)
  avals <- seq(first, last, length = k)

  for(i in 1:k) {
    xhat  <- diffeq((avals[i]-1), n = n, e = c(0, avals[i]*x))
    res   <- x - xhat
    SS[i] <- dot(res, res)
  }

  alpha <- avals[which(SS == min(SS))]
  SSmin <- SS[which(SS == min(SS))]

  xhat <- diffeq(alpha-1, n = n, e = c(0, alpha*x))
 
  del.x <- 0.05*(last - first)
  ylim <- c(min(x, xhat), max(x, xhat))
  par(mfrow = c(2, 1),  oma = rep(0, length = 4), mar = c(5, 4, 2, 2) + 0.1)

  plot(avals, SS, xlab = "a", ylab = "S(a)", pch = 20)
  mtext("Sum of Squares: Exponential Smoothing", line = 0.25, cex = 1.25)
  points(alpha, SSmin, col = "red", pch = 19)
  text(last-del.x, max(SS), labels = bquote(alpha == .(alpha)),  
       col = "red")
  plot(xhat, ylim = ylim, type = "p", pch = 20, col="blue", 
       xlab = "t", ylab = expression(x[t]))
  mtext("Series and Forecasts", line = 1.25, cex = 1.25)
  mtext(bquote(alpha == .(alpha)), line = 0.25, cex = 1)
  lines(1:n, x, lwd = 2)

  return(list(alpha = alpha, min.SS = SSmin, xhat = xhat, 
              avalues = avals, all.SS = SS)) 
}
#
#
#
extend <- function(x, h, d1, d2 = 0)
#-----------------------------------------------------------------------
#  R function to extend time series x by h points using antidifferencing.
#
#  INPUT: x = a real vector containing the series to be extended.
#         h = the number of points by which x should be extended.
#         d1 = the order of the first difference.
#         d2 = the order of the second difference.  If no second
#              differencing is needed, then d2 = 0 (the default).
#
#  VALUE: The function extend returns a vector of length n + h
#         containing the original series (as the first n elements) and
#         the extended series (in elements n+1, ..., n+h).
#-----------------------------------------------------------------------
{
  n <- length(x)
  if(d1 != 0 || h != 0) {
    x <- c(x, rep(0, length = h))
    if(d2 == 0) {
      zbar <- mean(diff(x[1:n], d1))
      for(i in 1:h) x[n+i] <- x[n+i-d1] + zbar
    } 
    else {
      zbar <- mean(diff(diff(x[1:n], d1), d2))
      for(i in 1:h) x[n+i] <- x[n+i-d1] + x[n+i-d2] - x[n+i-d1-d2] + zbar
    }	
  }
  return(x)
}
#
#
#
filt <- function(beta, beta0, x)
#-----------------------------------------------------------------------
#  R function to compute a filtered version of the vector x.
#
#  y_t = beta_0 x_{t+m} + beta_1 x_{t+m-1} + ... + \beta_m x_t,
#        t = 1, ..., n - m
#
#  INPUT: beta = a real vector containing the filter coefficients 
#                beta_1, \ldots, beta_m.
#         beta0 = a real scalar containing the filter leading
#                 coefficient.
#         x = a real vector to be filtered.
#
#  VALUE: The function filt returns a real vector of length n - m
#         containing the filtered version of x.
#
#  FORTRAN required: SUBROUTINE filt.
#-----------------------------------------------------------------------
{
  m <- length(beta)
  n <- length(x)
  y <- rep(0, length.out = n - m)
  if(m == 0) y <- beta0*x
  else {
    z <- .Fortran("filt",
                  as.double(beta),
                  as.double(beta0),
                  as.integer(m),
                  as.integer(n),
                  as.double(x),
              y = as.double(rep(0, length = n-m)))
  }
  return(z$y)
}
#
#
#
freqs <- function(n)
#-----------------------------------------------------------------------
#  R function to form a real vector of length [n/2] + 1 containing the 
#  natural frequencies between [0, 0.5] defined as  
#         freq[i] = (i - 1) / n,   i = 1 ,..., [n/2] + 1.
#
#  INPUT: n = an integer scalar containing the length of the original
#             time series.
#
#  VALUE: The function freqs returns a vector of length [n/2]+1 
#         containing the natural frequencies defined by
#         freq[i] = (i - 1) / n,   i = 1 ,..., [n/2] + 1.
#-----------------------------------------------------------------------
{
   return(c(0:(n/2))/n)
}
#
#
#
gram.schmidt <- function(x)
#-----------------------------------------------------------------------
#  R function to compute the Gram-Schmidt decomposition of a matrix x,
#  x = QR.
#
#  INPUT: x = a real matrix of dimension n x m.
#
#  VALUE: The function gram.schmidt returns a list containing the
#         following objects:
#         Q = real matrix of dimension n x m containing the Q factor in
#             the Gram-Schmidt decomposition.
#         R = real matrix of dimension m x m containing the R factor in 
#             the Gram-Schmidt decomposition.
#         error = an integer error indicator.  If error = 0, there was
#                 no error.  If error = 1, there was an error.
#
#  FORTRAN required: grmsmt SUBROUTINE
#-----------------------------------------------------------------------
{
  z <- .Fortran("grmsmt2ts",
                 as.double(x),
                 as.integer(nrow(x)),
                 as.integer(ncol(x)),
             Q = as.double(matrix(0, nrow(x), ncol(x))),
             R = as.double(matrix(0, ncol(x), ncol(x))),
         error = as.integer(0))

  error <- z$error
  if(error != 0) {
    Q <- NULL
    RR <- NULL
  }
  else {
    Q <- matrix(z$Q, nrow(x), ncol(x))
    RR <- matrix(z$R, ncol(x), ncol(x))
  }
  return(list(Q = Q, R = RR, error = error))
}
#
#
#
infqnt <- function(x)
#-----------------------------------------------------------------------
#  R function to construct the informative quantile plot of x.
#
#  INPUT: x = a real vector containing the data.
#
#  VALUE: The function infqnt creates an informative quantile plot of
#         the data in the real vector x.
#-----------------------------------------------------------------------
{
  n    <- length(x)
  xs   <- sort(x)
  infq <- clip((xs - median(xs))/(2*(xs[3*n/4] - xs[n/4])), -1, 1)
  vs <- ((1:n) - 0.5) / n

  plot(vs, infq, type = "l", lwd = 3, ylim = c(-1, 1), 
       xlab = expression(u[(i)]), ylab = expression(y[(i)]),
       main = "Informative Quantile Plot", frame = FALSE)
  lines(c(0, 1), c(-0.5, 0.5))
  lines(c(0, 0, 1, 1, 0), c(-1, 1, 1, -1, -1))
#
  return(invisible())
}
#
#
#
macorr <- function(beta, m = length(beta), rvar = 1)
#-----------------------------------------------------------------------
#  R function to evaluate autocorrelation function and variance for an  
#  MA process with MA coefficients beta and error variance rvar.
#
#  INPUT: beta = a real vector of length q containing the MA coefficients.
#         m = an integer scalar containing maximum autocorrelation lag.
#             The default value is m = length(beta).
#         rvar = a real scalar containing error variance.  The default
#                is rvar = 1.
#
#  VALUE: The function arcorr returns a list containing the following
#         objects:
#         var = a real scalar containing the variance of the MA 
#               process.
#         corr = a real vector of length m containing the
#                autocorrelation function of the MA process.
#         error = an integer scalar that is an error indicator.  If 
#                 error = 0, there was no error.  If error = 1, then
#                 the value of p or q was negative.  If error = 2,
#                 then the value of m was less than max(p,q) + 1.
#                 If error = 3, then the MA model was not stationary.
#
#  FORTRAN required: SUBROUTINE wilson
#-----------------------------------------------------------------------
{
  return(armacorr(alpha = 0, beta, m, rvar))
}
#
#
#
madt <- function(beta, n, rvar = 1, seed = 0)
#----------------------------------------------------------------------
#  R function to simulate a realization of length n from an MA process
#  having coefficients beta and error variance rvar.
#
#  INPUT: beta = a real vector containing the MA coefficients.
#         n = an integer scalar specifying the length of the desired
#             realization.
#         rvar = an real scalar containing the variance of the errors.
#                The default value is rvar = 1.
#         seed = a nonnegative integer for the random number generator.
#                If seed = 0 (the default) then the output seed from
#                the last call is used.  "seed" is an optional argument.
#
#  VALUE: The function madt returns a vector of length n containing 
#         the desired realization.
#----------------------------------------------------------------------
{
   eps <- sqrt(rvar)*wn(n = (n + length(beta)), seed = seed)

   return(filt(beta = beta, beta0 = 1, x = eps))
}
#
#
#
mapart <- function(beta, rvar = 1, m = length(beta) + 10)
#-----------------------------------------------------------------------
#  R function to compute MA partial partial autocorreltions.
#
#  INPUT: beta = a real vector of length q containing the MA 
#                coefficients.
#         rvar = a real scalar containing the error variance.  The
#                default value is rvar = 1.
#         m = an integer scalar specifying the number of partial
#             autocorrelations to return.  The default value is 
#             m = q + 10.
#
#  VALUE: The function armapart returns a list containing the following
#         objects:
#         theta = a real vector of length m containing the partial
#                 autocorrelations of the MA process.  If error is
#                 not zero, then theta is an empty vector.
#         error = an integer scalar indicating if there was no error
#                 (error = 0) or if the specified model was 
#                 nonstationary (error not 0).
#------------------------------------------------------------------------
{
  return(armapart(alpha = 0, beta, rvar, m))
}
#
#
#
masmooth <- function(x, k)
#-----------------------------------------------------------------------
#  R function to apply a moving average smoother of length k to x.
#
#  INPUT: x = a real vector containing the series to be smoothed.
#         k = an integer scalar containing the order of the moving
#             average smoother.
#
#  VALUE: The function masmooth returns a vector of containing the
#         smoothed version of the series x.
#-----------------------------------------------------------------------
{
  k <- as.integer(k)

  z <- x
  if(k != 0) {
    n <- length(x)
    w <- c(x[(k+1):2], x, x[(n-1):(n-k)])
    for(i in 1:k)  z <- z + w[(k+1-i):(k-i+n)] + w[(k+1+i):(k+i+n)]
    z <- z / (2 * k + 1)
  }
  return(z)
}
#
#
#
masp <- function(beta, rvar = 1, Q = 256)
#-----------------------------------------------------------------------
#  R function to compute the spectral density function at Q frequencies
#  in [0, 1] for an MA(q) process with coefficients beta and error 
#  variance rvar.
#
#  INPUT: beta = a real vector of length q containing the coefficients
#                of the MA process.
#         rvar = a real scalar containing the error variance.  The 
#                default value is rvar = 1.
#         Q = an integer scalar indicating the number of frequencies
#             in the interval [0, 1] at which to compute the spectral 
#             density function.  The default value is Q = 256.
#
#  VALUE: The function masp returns a real vector of length [Q/2]+1
#         containing the spectral density function of the MA(q) 
#         process.
#-----------------------------------------------------------------------
{
  return(armasp(alpha = 0, beta, rvar, Q))
}
#
#
mchol <- function(x)
#-----------------------------------------------------------------------
#  R function to compute the modified Cholesky decomposition of a real
#  symmetrix matrix x: find lower triangular matrix L and diagonal 
#  matrix D such that x = LDL'.
#
#  INPUT: x = a real square (dimension m) symmetrix matrix.
#
#  VALUE: The function mchol returns a list containing the following
#         objects:
#         L = a real m x m lower triangular matrix in the modified
#             Cholesky decomposition.
#         D = a real m x m diagonal matrix in the modified Cholesky
#             decomposition.
#         error = an error indicator.  If error = 0 there was no error.
#                 If error = 1, the modified Cholesky decomposition 
#                 could not be found.
#
#  FORTRAN required: SUBROUTINE mchol2ts
#-----------------------------------------------------------------------
{
  if(nrow(x) != ncol(x)) stop("\n x is not a square matrix.\n")
  if(max(y - t(y)) > 1e-10) stop("\n x is not a symmetrix matrix.\n")

  m <- nrow(x)
  error <- 0

  z <- .Fortran("mchol2ts", 
                as.double(x),
                as.integer(m), 
            L = as.double(matrix(0, nrow = m, ncol = m)),
            D = as.double(rep(0, length = m)),
        error = as.integer(error))

  error <- z$error
  if(error == 1) {
    L <- NULL
    D <- NULL
  }
  else {
    L <- t(matrix(z$L, nrow(x), nrow(x)))
    D <- diag(z$D)
  }
  return(list(L = L, D = D, error = z$error))
}
#
#
#
multpoly <- function(alpha, beta)
#-----------------------------------------------------------------------
#  R function to multiply two polynomials.
#
#  INPUT: alpha = a real vector containing the coefficients of one of
#                 the two polynomials.
#         beta = a real vector containing the coefficients of the
#                second polynomial.
#  NOTE: The coefficient of the lead-order term are assumed to be one
#        for both polynomials and should not be included in alpha or
#        beta.
#
#  VALUE: The function multpoly returns a real vector containing the 
#         coefficients of the product of the two polynomials, beginning
#         with the (p+q)-1 order term (since the coefficient of the
#         p+q-order term is one).
#
#  FORTRAN required: SUBROUTINE mtpoly
#
#  EXAMPLE: To multiple x^2 + 3 x - 1 and x - 3, the arguments are
#           alpha = c(3,-1), beta = -3.
#           The function returns 0  -10  3.  So the product is
#           x^3 - 10 x + 3.
#-----------------------------------------------------------------------
{
  p <- length(alpha)
  q <- length(beta)

  if(p == 0 || q == 0) 
    stop("\nThe order of both polynomials must be greater than zero.\n")

  z <- .Fortran("mtpoly",
                as.double(alpha),
                as.double(beta),
                as.integer(p),
                as.integer(q),
        gamma = as.double(rep(0, p+q)))

  return(z$gamma)
}
#
#
#
parcorr <- function(x, m) 
#-----------------------------------------------------------------------
#  R function to compute the sample partial autocorrelation function
#  function of a time series (x) for lags 1, ..., m.
#
#  INPUT: x = a real vector containing the time series.
#         m = an integer scalar indicating the number of lags for
#             computing the sample partial autocorrelation function.
#
#  VALUE: The function parcorr returns a real vector of length m
#         containing the sample partial autocorrelation function of x.
#-----------------------------------------------------------------------
{
   e     <- c( x - mean(x), rep(0, length = m) )
   f     <- crlag(e)
   theta <- rep(0, length = m)

   for(i in 1:m) {
      part     <- dot(f, e) / dot(f, f)
      temp     <- e
      e        <- temp - part * f
      f        <- crlag(f - part * temp)
      theta[i] <- part
   }
#
   return(theta)
}
#
#
#
parz <- function(m)
#-----------------------------------------------------------------------
#  R function to compute the Parzen weight function for indices 
#  0, 1, ..., m (on frequencies 0, 1/m,..., (m-1)/m, 1).
#
#  INPUT: m = an integer scalar specifying the number of indices at 
#             which to compute the Parzen weights.
#
#  VALUE: The parz function returns a real vector of length m 
#         containing the value of the Parzen weight function for indices 
#         0, 1, ..., m (on frequencies 0, 1/m, ..., (m-1)/m, 1).
#-----------------------------------------------------------------------
{
  u  <- (0:m) / m
  z1 <- u[which(u <= 0.5)]
  z2 <- u[which(u > 0.5)]
  z1 <- 1 - 6*z1*(z1 - z1*z1)
  z2 <- 2*(1 - z2)^3

  return(c(z1,z2))
}
#
#
#
perdave <- function(n, alpha, nsamps = 20, Q = n, rvar = 1, seed = 0)
#-----------------------------------------------------------------------
#  R function to illustrate the unbiasedness of the periodogram.  For
#  a specified AR(p) process, the function generates nsamps samples of
#  size n.
#
#  INPUT: n = an integer scalar indicating the length of the each
#             series.
#         alpha = a real vector containing the coefficients of the AR
#                 process.
#         nsamps = the number of samples to generate.  The default
#                  value is nsamps = 20.
#         Q = an integer specifying the number of frequencies on the
#             interval [0, 1] at which to compute the spectral density
#             function.  The default value is Q = n (which is best for
#             comparing returned numerical output).
#         rvar = a real scalar containing the variance of the error
#                terms.  The default value is rvar = 1.
#         seed = a nonnegative integer for the random number generator.
#                If seed = 0 (the default) then the output seed from
#                the last call is used.  "seed" is an optional argument.
#
#  VALUE: The function perdave creates a plotting window with four
#          graphs in a 2 x 2 array.  The top left graph is the true
#          AR(p) spectral density function.  The top right graph is
#          the average of the periodograms generated thus far.  The
#          bottom left graph is of the current AR(p) realization.
#          The bottom right graph is of the periodogram of the
#          current AR process.  When the function exits, a list 
#          containing the following objects is returned.
#          truef = a vector of length [Q/2]+1 containing the true
#                  AR spectra,
#          averagef = a vector of length [n/2]+1 containing the
#                     average of the periodograms.
#          nsamps = the number of samples used for computing averagef.
#-----------------------------------------------------------------------
{
  p        <- length(alpha)
  R0       <- arcorr(alpha, rvar, m = (p + 1))$var
  f        <- arsp(alpha, rvar, Q)
  averagef <- rep(0, length = (floor(n/2)+1))

  par(mar = c(3, 2, 2, 1), oma = rep(0, length = 4),
      tck = -0.015, cex.axis = 0.75, xaxs = "i", yaxs = "i",
      cex.main = 0.75, mgp = c(1, 0.25, 0), cex.lab = 0.85)

  split.screen(c(2, 2), erase = TRUE)
  screen(n = 1)
  plotsp(f, Q, R0, main = paste("True AR(",p,") Spectral Density"))

  samp <- 0
  b    <- "C"
  while(b == "C" && samp < nsamps) {
    samp     <- samp + 1
    x        <- ardt(alpha, n, rvar, seed)$x
    fhat     <- perdgm(x)
    averagef <- ((samp - 1)*averagef + fhat) / samp

    screen(n = 2, new = TRUE)
    plotsp(averagef, n, R0, main = paste("Average of",samp,"Periodograms"))
    screen(n = 3, new = TRUE)
    plot(x, type = "l", xlab = "t", ylab = expression(x[t]),
         main = paste("AR(",p,") Realization #",samp))
    screen(n = 4, new = TRUE)
    plotsp(fhat, n, R0, main = paste("Periodogram of Realization #",samp))

    if(samp < nsamps) {
      b <- readline("Do you want to continue? (Strike <Enter>, or type C or Y.) ")
      if(b == "c" || b == "Y" || b == "y" || b == "") b <- "C"
    }
  }
  if(samp < nsamps) nsamps <- samp
  close.screen(all = TRUE)

  return(list(f = f, averagef = averagef, nsamps = nsamps))
}
#
#
#
perdgm <- function(x)
#-----------------------------------------------------------------------
#  R function to compute the periodogram of the time series x.
#  Priestleym p. 395
#
#  INPUT: x = a real vector containing the time series (of length n).
# 
#  VALUE: The function perdgm returns a vector of length [n/2]+1
#         containing the periodogram t the natural frequencies.
#-----------------------------------------------------------------------
{
   ft <- 2*Mod(fft(x - mean(x), inverse = TRUE))^2 / (n <- length(x))
   ft <- ft[1:((n / 2) + 1)]
   return(ft)
}
#
#
#
phaseplot <- function(x, k, main = "Phase Plot")
#-----------------------------------------------------------------------
#  R function to create a scatterplot of x_{t-k} versus x_t.
#
#  INPUT: x = a real vector containing the time series.
#         k = an integer scalar containing the lag.
#         main = a character string containing the title for
#                the plot.
#
#  VALUE: The function phaseplot creates a scatterplot of
#         x_{t-k} verssu x_t.
#-----------------------------------------------------------------------
{
  plot(x[1:(length(x)-k)], x[(k+1):length(x)], pch = 20, main = main,
       xlab = expression(x[t]), ylab = bquote(x[t-.(k)]))

  return(invisible())
}
#
#
#
plotcorr <- function(r, m = length(r), corr.type = "acf", type = "l",
                     main = "Autocorrelation Function")
#----------------------------------------------------------------------
#  R function to create a plot of the autocorrelation function (acf) or
#  partial autocorrelation function (pacf).
#
#  INPUT: r = a real vector containing the values of the autocorre-
#             lation or partial autocorrelation function.
#         m = an integer scalar indicating how many lags to plot the
#             acf or pacf.  The default is m = length(r).  If 
#             m > length(r), then the length(r)+1, ..., m values of
#             r will be assumed zero.
#         corr.type = a character string indicating whether the values
#                are autocorrelations or partial autocorrelations.  The
#                default is type = "acf".  For partial autocorrelation,
#                use type = "pacf".
#         type = a character indicating what type of plot should be
#                drawn.  The default is type = "l". Possible types are
#                "p" for points, "l" for lines, "b" for both points
#                and lines, or "v" for vertical lines with points.
#         main = a character string containing the title of the plot.
#                The default is main = "Autocorrelation Function".
#
#  RESULT: The function plotcorr creates a plot of the autocorrelation
#          or partial autocorrelations at lags 1, 2, ..., m contained
#          in r.  No numerical output is returned.
#
#  Written: 11/18/2015, JLH
#----------------------------------------------------------------------
{
  if(corr.type == "acf")  ylab <- expression(rho[nu])
  if(corr.type == "pacf") ylab <- expression(theta[nu])
  
  if(m > length(r)) {
    r <- c(r, rep(0, times = (m - length(r)))) 
  }
  else {
    r <- r[1:m]
  }
  
  plot(1:m, r, ylim = c(-1, 1), type = "n", xlab = expression(nu),
       ylab = ylab, main = main)
  abline(h = 0)
  if(type == "p" || type == "v") points(r, pch = 19)
  if(type == "v") {
    for(j in 1:m) lines(c(j, j), c(0, r[j]))
  }
  if(type == "l") lines(r)
  
  return(invisible())
}
#
#
#
plotsp <- function(f, n, div, 
                   main = "Natural Log of Standardized Spectrum")
#----------------------------------------------------------------------
#  R function to plot the natural logarithm of the standardized version
#  of the estimated spectral density function in f.
#
#  INPUT: f = a real vector of length [n/2]+1 containing the estimated
#             spectral density function.
#         n = an integer scalar containing the length of the original
#             series.
#         div = a real scalar containing the divisor for standardizing
#               the estimated spectral density (in f).
#         main = a character string containing the title for the
#                plot.
#
#  VALUE: The function plotsp creates a plot of the natural logarithm
#         of the standardized estimated spectral density function.
#-----------------------------------------------------------------------
{
  plot(freqs(n), log(stdf(f, div, exp(-6), exp(6))), type = "l",
       xlab = expression(omega), ylab = expression(hat(f)*(omega)),
       main = main, xlim = c(0, 0.5), ylim = c(-6, 6))
#
  return(invisible())
}
#
#
#
polar <- function(z = NULL, a = NULL, b = NULL) 
#-----------------------------------------------------------------------
#  R function to find the amplitude (modulus) and phase of a complex
#  number z.  If z is not provided, then a is the real part of the
#  complex number and b is the imaginary part of the complex number.
#
#  INPUT: Either z = a complex number (vector or scalar) OR
#         BOTH   a = the real part of the complex number and
#                b = the imaginary part of the complex number.
#  If z is provided, a and b are ignored.
#
#  VALUE: The function polar returns a list containing the the following
#         two objects with the same dimension as z (or as a and b).
#         A = a real array containing the the amplitude of the complex 
#             number(s), and
#         phi = a real array containg the phase of the complex number(s).
#-----------------------------------------------------------------------
{
  if(is.null(z) && is.null(a) && is.null(b))
    stop("\nProvide either z or a and b.\n")
  if(is.null(z) && length(a) != length(b))
    stop("\nThe lengths of a and b must be equal.\n")

  if(is.null(z)) z <- complex(real = a, imaginary = b)
  if(!is.complex(z)) z <- as.complex(z)
  if(is.complex(z)) {
    a <- Re(z); b <- Im(z)
  }
#  
  A <- sqrt(a*a + b*b)
#
  phi <- rep(0, length = length(a))
  inds <- which(a < 0)
  phi <- atan(b/a)
  phi[inds] <- phi[inds] + pi
 
  return(list(A = A, phi = phi))
}
#
#
#
poly <- function(coeffs, x)
#-----------------------------------------------------------------------
#  R function to evaluate a polynomial having coefficients coeffs 
#  arranged from degree 0, 1, ..., p at the values in the vector x.
#
#  INPUT: coeffs = a real vector containing the coefficients of the
#                  polynomials.
#         x = a real vector containing the values at which to evaluate
#             the polynomial.
#
#  VALUE: The function poly returns a real vector containing the
#         values at x of the polynomial with coefficients coeffs.
#
#  FORTRAN required: SUBROUTINE poly
#-----------------------------------------------------------------------
{
  deg <- length(coeffs) - 1
  n   <- length(x)

  if(deg == 0) 
    fx <- rep(coeffs[1], length(x))
  else {
    z <- .Fortran("poly",
                  as.double(coeffs),
                  as.integer(deg),
                  as.double(x),
                  as.integer(n),
             fx = as.double(rep(0, length = n)))
    fx <- z$fx
  }
  return(fx)
}
#
#
#
polyrt <- function(coeffs, m = 100, eps = 1.e-6)
#-----------------------------------------------------------------------
#  R function to compute the roots (zeros) of a polynomials
#  1 + sum_{j=1}^p coeffs_j z^j given its coefficients.
#
#  INPUT: coeffs = a vector containing the coefficients of the
#                  polynomial.  It is assumed the constant terms is one.
#         m = an integer scalar containing the number of iterations.
#         eps = a real scalar containing a convergence criterion.
#
#  VALUE: The function polyrt returns a list containing the following
#         objects:
#         roots = a possibly complex vector of length p containing
#                 the zeros (roots) of the polynomial.
#         error = an integer error indicator.  If error = 0, the
#                 algorithm converged within eps in no more than m
#                 iterations.  Otherwise, error > 0.
#
#  FORTRAN required: SUBROUTINE polyrt
#-----------------------------------------------------------------------
{
  p     <- length(coeffs)
  error <- 0

  z <- .Fortran("polyrt",
                as.double(coeffs),
                as.integer(p),
                as.integer(m),
                as.double(eps),
        roots = as.double(matrix(0,2,p)),
        error = as.integer(error))
  
  error <- z$error
  if(error != 0) roots <- NULL
  if(error == 0) {
    roots <- matrix(z$roots, 2, p)
    roots <- complex(real = roots[1, 1:p], imag = roots[2, 1:p])
  }
  return(list(roots = roots, error = error))
}
#
#
#
rtpoly <- function(roots)
#-----------------------------------------------------------------------
#  R function to find the coefficients 1 + sum_{j=1}^p a_j z^j of a
#  polynomial, given its roots.
#
#  INPUT: roots = a possibly complex vector containing the roots of the 
#                 polynomial.
#
#  VALUE: The function rtpoly returns a vector of length p containing
#         the coefficients a_1, ... a_p of a polynomial having roots
#         given in roots.
#
#  FORTRAN required: SUBROUTINE rtpoly
#-----------------------------------------------------------------------
{
  p <- length(roots)
  if(!is.complex(roots)) roots <- as.complex(roots)

  z <- .Fortran("rtpoly",
                as.double(rbind(Re(roots), Im(roots))),
                as.integer(p),
       coeffs = as.double(rep(0, length = p)))

  return(z$coeffs)
}
#
#
#
rw <- function(n, seed = 0) 
#-----------------------------------------------------------------------
#  R function to generate a random walk of length n using normal
#  white noise.
#
#  INPUT: n = an integer scalar containing the length of the desired
#             random walk.
#         seed = a nonnegative integer for the random number generator.
#                If seed = 0 (the default) then the output seed from
#                the last call is used.  "seed" is an optional argument.
# 
#  VALUE: The function rw returns a realization of length n from a
#         Gaussian random walk.
#-----------------------------------------------------------------------
{
   if(seed != 0) set.seed(seed)
   return(cumsum(rnorm(n)))
}
#
#
#
srvars <- function(x, m) 
#-----------------------------------------------------------------------
#  R function to compute m standardized residual variances of a time
#  series x.
#
#  INPUT: x = a real vector containing the time series.
#         m = an integer scalar containin the number of lags for
#              computing the standardized residual variances.
#
#  VALUE: The function srvars returns a real vector of length m 
#         containing the standardized residual variances of x.
#------------------------------------------------------------------------
{
  if(m == 0) xx <- 0 else {
    part <- parcorr(x, m)
    xx   <- exp(cumsum(log(1 - part*part)))
  }
#
  return(xx)
}
#
#
#
submns <- function(x, d = 1, option = 1, xbar = rep(0, length.out = d))
#----------------------------------------------------------------------
#  R function to subtract out seasonal means from a time series x.
#
#  INPUT: x = a real vector containing the time series.  The length
#             of x must be a multiple of d.
#         d = an integer scalar indicating the length of the cycle;
#             for example, for monthly data across a year, d = 12.
#             The default value of d = 1.  If d = 1, the function
#             computes the sample mean of the series.
#         option = an integer scalar indicating what operation should
#             be performed.  The default value of option is option = 1.
#         xbar = a vector of length d containing the seasonal means.
#                The default value of xbar is rep(0, times = d). 
#                xbar is ignored if option = 0 or 1.
#
#  The result of submns is dependent upon the value of the option
#  argument.
#  If option = 0: the function will return a list containing 
#       x = a real vector containing the original series, and
#       xbar = a vector of length d containing the seasonal means.
#     option = 1: the function will return a list containing
#       x = a real vector containing the original series with the
#             seasonal means subtracted. 
#       xbar = a vector of length d containing the seasonal means.
#     option = 2: the function will return the series x with the
#       seasonal means provided in xbar added to the series.
#----------------------------------------------------------------------
{
  n <- length(x)
  y <- x
  if(d == 1) {
    if(option == 0) return(mean(x))
    if(option == 1) return(list(xbar= mean(x),x = x-mean(x)))
    if(option ==2 ) return(x+xbar)}
  
  if(option == 2) {
    for(i in 1:d) {
      inds <- seq(from = i, to = n, by = d)
      x[inds] <- x[inds] + xbar[i]}
    return(x)
  }
  
  xbar <- rep(0, times = d)
  for(i in 1:d) {
    inds <- seq(from = i, to = n, by = d)
    xbar[i] <- mean(x[inds])
    if(option == 1) x[inds] <- x[inds] - xbar[i]
  }
  return(list(x = x, xbar = xbar))
}
# 
#
#
stdf <- function(f, fac, a, b) 
#-----------------------------------------------------------------------
#  R function to divide f by fac and clip the results to be between
#  the real scalars a and b.
#
#  INPUT: f = a real vector.
#         fac = a real scalar containing the dividing factor.
#         a = a real scalar containing the lower bound for the clip.
#         b = a real scalar containing the upper bound for the clip.
#
#  VALUE: The function stdf returns a real vector containing the
#         values f/fac clipped so that all values are between a and b.
#-----------------------------------------------------------------------
{
   f   <- f / fac
   f1  <- c(1:length(f))
   f   <- replace(f, f1[f > b], b)
   f   <- replace(f, f1[f < a], a)
#
   return(f)
}
#
#
#
swp <- function(A, k1 = 1, k2 = ncol(A))
#---------------------------------------------------------------------
#  R function to sweep the square matrix A on diagonals k1, ..., k2.
#
#  INPUT: A = a real square matrix to be swept.
#         k1 = an integer containing the index of the first diagonal 
#              of A on which to begin sweeping.  The default value is
#              k = 1.
#         k2 = an integer containing the index of the last diagonal 
#              of A on which to end sweeping.  The default value is
#              k = ncol(A).
#
#  VALUE: The function swp returns a list containing the following
#         objects:
#         A = the result of sweeping the input matrix on diagonals
#             k1 through k2.
#         error = an integer error indicator.  If error = 0, then no
#                 error occurred.  Otherwise, error = 1.
#
#  FORTRAN required: SUBROUTINE swpk12
#----------------------------------------------------------------------
{
  if(ncol(A) != nrow(A)) stop("\n The matrix A must be square.\n")
  if(k1 < 1 || k2 > n || k1 > k2) stop("\nIllegal k1 or k2 in swp().\n")

  n <- ncol(A)

  z <- .Fortran("swpk12", 
           a1 = as.double(A), 
                as.integer(n),
                as.integer(k1),
                as.integer(k2), 
          ier = as.integer(0))

  error <- z$ier
  if(error == 0) a1 <- matrix(z$a1, n, n)
  else a1 <- NULL

  return(list(A = a1, error = z$ier))
}
#
#
#
toepl <- function(R, R0, m = (length(R)+1))
#-----------------------------------------------------------------------
#  R function to form a Toeplitz matrix, Toepl(R0, ..., R(m-1)).
#
#  INPUT: R = a real vector containing off-diagonal elements of the
#             desired Toeplitz matrix.
#         R0 = a real scalar containing the diagonal element of the
#              desired Toeplitz matrix.
#         m = an integer scalar containing the dimension of the desired
#             Toeplitz matrix.  The default is length(R) + 1.  If a 
#             lesser value of m is specified, then the function will use
#             the first m-1 elements of R.
#
#  VALUE: The function toepl returns a Toeplitz matrix of dimension m
#         with diagonal elements R0 and off-diagonal elements in R.
#-----------------------------------------------------------------------
{
  if(m > length(R)+1) stop("\nm must be less than one plus the length of R.\n")
  G <- matrix(0, m, m)
  r <- c(R0, R[1:(m-1)])

  return(matrix(r[abs(col(G) - row(G)) + 1], m, m))
} 
#
#
#
tsdensity <- function(x, rbins = 0, npts = 501, kernel = "parzen") 
#-----------------------------------------------------------------------
#  Timeslab in R function to compute kernel density estimate of 
#  univariate data.
#
#  INPUT: x = a real vector containing the data
#         rbins = a real scalar indicating the method for computing
#                 the bandwidth.  If rbins > 0, then h = r_x/(2*rbins).
#                 Otherwise, rbins = r_x n^q/(3.5 s_x). q = 1/3 for
#                 rectangular window and 1/5 for all others.
#         npts = an integer scalar indicating the number of points at
#                which to compute the kernel density estimate.
#         kernel = a character or integer argument indicating which
#                  kernel function to use.  The default is
#                  kernel = "parzen" (or kernel = 7).  Other values are
#                  and the associated kernels are
#                  "rectangular" or 1       "biweight" or 4
#                  "cosine" or 2            "triangular" or 5
#                  "epanechnikov" or 3      "gaussian" or 6
#                  The argument may be truncated to the first (unique)
#                  character of the name of the window.
#
#  VALUE: The function tsdensity returns a list containing the following
#         objects.
#         kernel = the input for the argument "kernel",
#         bw = the bandwidth used in smoothing,
#         rbins = the input for the argument "rbins", or the computed
#                 value of rbins if the input was <= 0,
#         y = a real vector of length npts contains the set of 
#             frequencies (x's) at which the density is estimated, and 
#         fy = a real vector of length npts containing the kernel
#              density estimate for the corresponding element in y.
#
#  Written: 10/12/2010 JLH
#-----------------------------------------------------------------------
{
  kopt  <- 0
  kfunc <- substring(as.character(kernel), 1, 1)
  if(kfunc == "r" || kfunc == "R" || kfunc == "1") kopt <- 1
  if(kfunc == "c" || kfunc == "C" || kfunc == "2") kopt <- 2
  if(kfunc == "e" || kfunc == "E" || kfunc == "3") kopt <- 3
  if(kfunc == "b" || kfunc == "B" || kfunc == "4") kopt <- 4
  if(kfunc == "t" || kfunc == "T" || kfunc == "5") kopt <- 5
  if(kfunc == "g" || kfunc == "G" || kfunc == "6") kopt <- 6
  if(kfunc == "p" || kfunc == "P" || kfunc == "7") kopt <- 7

  if(kopt == 0) stop("\nIllegal input for argument kernel.\n")

  if(npts < 0)  npts <- -npts
  if(npts == 0) npts <- 501

  rx <- diff(range(x))
  if(rx < 1e-10) stop("\n The range of the data is too small.\n")

  if(kopt == 1) q <- 1/3 else q <- 1/5
  if(rbins <= 0) rbins <- rx*npts^q/(3.5*sd(x)) + 1
  h <- rx/(2*rbins)
 
  y    <- seq(min(x), max(x),  length = npts)
  fhat <- rep(0, length = npts)
  n    <- length(x)

  fy <- .Fortran("tsdensity",
                 as.integer(n),
                 as.integer(npts),
                 as.integer(kopt),
                 as.double(h),
                 as.double(x),
                 as.double(y),
          fhat = as.double(fhat))

  return(list(bw = h, y = y, fy = fy$fhat))
}
tsplot <- function(x, xlab = NULL, ylab = NULL, main = "Time Series Plot") 
#----------------------------------------------------------------------
# R function to plot a time series.
#
# INPUT: x = a real vector containing the series to be plotted.
#        xlab = a character string containing the label for the time
#               (horiztonal) axis. If no value is supplied, then 
#               xlab = "t".
#        ylab = a character string containing the label for the 
#               vertical axis.  If no value is supplied, then
#               ylab = expression(x[t]).
#        main = a character string containing the title for the plot.
#               If no value is supplied, then main = "Time Series Plot".
#
# OUTPUT: The function tsplot creates a time series plot of the series x,
#         having time axis label xlab and vertical axis label ylab, with
#         title main.  No numerical output is returned.
#----------------------------------------------------------------------
{
  if(is.null(xlab)) xlab <- "t"
  if(is.null(ylab)) ylab <- expression(x[t])
  plot(x, type = "l", xlab = xlab, ylab = ylab, main = main)
  return(invisible())
}
#
#
#
tsvar <- function(x) 
#-----------------------------------------------------------------------
#  R function to find the sample variance of a time series.
#
#  INPUT: x = a real vector containing the time series.
#
#  VALUE: The function tsvar returns a real scalar equal to the sample
#         variance of the time series.
#-----------------------------------------------------------------------
{
   return(dot(x - mean(x), x - mean(x)) / length(x))
}
#
#
#
wn <- function(n, dist = "normal", seed = 0) 
#-----------------------------------------------------------------------
#  R function to simulate white noise.
#
#  Generate a white noise time series of length n having distribution
#  specified by dist as described in Table 1.1 "Distributions of the wn 
#  function" of the Timeslab text.
#
#  INPUT:  n = an integer scalar containing the desired length of the
#              white noise series.
#          dist = an integer or character valued argument indicating
#                 the distribution of the white noise.  Values of dist
#                 are supplied in the Timeslab text.  If the name of
#                 the distribution is used, it must be as specified in
#                 the table (in the text).  The default is 
#                 dist = "normal".  Other permissible values are
#                 dist = 1 or "normal" ==> normal 
#                 dist = 2 or "uniform" ==> uniform 
#                 dist = 3 or "exponential" ==> exponential 
#                 dist = 4 or "logistic" ==> logistic 
#                 dist = 5 or "cauchy" ==> Caucy 
#                 dist = 6 or "extvalue" ==> extreme value
#                 dist = 7 or "lognormal" ==> lognormal
#                 dist = 8 or "dexponential" ==> double exponential
#                 If the character value is used, then the input may
#                 truncated to the first FOUR unique characters of
#                 the full string.
#          seed = a nonnegative integer for the random number generator.
#                 If seed = 0 (the default) then the output seed from
#                 the last call is used.  "seed" is an optional argument.
#
#  VALUE: The function wn returns a vector of length n containing
#         a realization of white noise from the specified distribution.
#------------------------------------------------------------------------
{
  if(seed[1] != 0) set.seed(seed)

  dfunc <- substring(as.character(dist), 1, 4)
  dopt <- 0
  if(dfunc == "norm" || dfunc == "Norm" || dfunc == "1") dopt <- 1
  if(dfunc == "unif" || dfunc == "Unif" || dfunc == "2") dopt <- 2
  if(dfunc == "expo" || dfunc == "Expo" || dfunc == "3") dopt <- 3
  if(dfunc == "logi" || dfunc == "Logi" || dfunc == "4") dopt <- 4
  if(dfunc == "cauc" || dfunc == "Cauc" || dfunc == "5") dopt <- 5
  if(dfunc == "extv" || dfunc == "Extv" || dfunc == "6") dopt <- 6
  if(dfunc == "logn" || dfunc == "Logn" || dfunc == "7") dopt <- 7
  if(dfunc == "dexp" || dfunc == "Dexp" || dfunc == "7") dopt <- 8

  if(dopt == 0) stop("\nIllegal input for argument dist.\n")

  if(dopt == 1) u <- rnorm(n)

  if(dopt == 2) u <- runif(n)

  if(dopt == 3) u <- -log(1 - u)

  if(dopt == 4) u <- log(u/(1 - u))

  if(dopt == 5) u <- tan(pi*(u - 0.5))

  if(dopt == 6) u <- log(-log(1 - u))

  if(dopt == 7) u <- exp(rnorm(n))

  if(dopt == 8) { 
    xx    <- (1:n)[u > 0.5]
    u[xx] <- 1 - u[xx] 
    u     <- log(2*u)
    u[xx] <- -u[xx]
  }
  return(u)
}
#
#
#
wntest <- function(x, m = floor(10*log10(length(x))), conf = 0.95)
#-----------------------------------------------------------------------
#  R function to find approximate (1-conf)100% simultaneous 
#  confidence bands for the correlogram and cumulative periodogram
#  under the null hypothesis of white noise.
#
#  Input: x = a real vector containing the series.
#         m = an integer scalar containing the maximum lag for the
#             sample autocorrelation function.  The default is 
#             m = [10 log_{10}(n)], where n = length(x).
#         conf = a real scalar containing the desired level of
#                confidence for the simultaneous confidence bands.
#                The default is conf = 0.95.  The value of conf must
#                be between 0.5 and 1.
#
#  VALUE: The function wntest creates a set of two graphs, in a 2x1 
#         array. The top graph contains the correlogram with the 
#         simultaneous confidence bands superimposed.  The bottom graph 
#         contains the cumulative periodogram with simultaneous 
#         confidence bands superimposed.
#-----------------------------------------------------------------------
{  
   if(conf < 0.5 || conf > 1) 
     stop("\nThe value of conf must be between 0.5 and 1.\n")
   if(m >= length(x) || m < 1)
     stop("\nThe value of m must be greater than 0 and less than n-1.\n")

   simul.alpha <- (1 + conf^(1/m))/2
   rho.hat     <- corr(x, m)$corr
   cum.perd    <- cumsp(perdgm(x))
   upper.rho   <- qnorm(simul.alpha)/sqrt(length(x))
   lower.rho   <- -upper.rho
   simul.alpha <- round(simul.alpha, 6)
#
   fb <- rep(0, 501)
   b  <- seq(0.4, 2, length = 501)
   for(i in 1:501) fb[i] <- 1 - barttest(b[i], length(x))$pval
   inds <- which(round(fb, 3) == round(conf, 3))
   if(length(inds) == 1) cc <- b[inds] else cc <- median(b[inds])
   xx <- seq(0, 0.5, length = length(x))
   yy.upper <- 2*xx + cc/sqrt(length(cum.perd))
   yy.lower <- 2*xx - cc/sqrt(length(cum.perd))
#
   par(mfrow = c(2, 1), oma = rep(0, length = 4), mar = c(5, 4, 4, 2) + 0.1)
   plot(1:m, rho.hat, type = "l", ylim = c(-1, 1), xlab = "v",
        ylab = expression(hat(rho[v])), 
        main = "Correlogram for White Noise Test")
   lines(c(1, m), rep(upper.rho, 2), lty = 2)
   lines(c(1, m), rep(lower.rho, 2), lty = 2)
   plot(freqs(length(x)), cum.perd, ylim = c(0, 1), pch = 20,
        xlab = expression(omega), ylab = expression(hat(F)*(omega)),
        main = "Cumulative Periodogram Test for White Noise")
   lines(xx, yy.upper, lty = 2)
   lines(xx, yy.lower, lty = 2)
#
   return(invisible())
}
#
#
#
datasets <- function()
#-----------------------------------------------------------------------
#  R function to assign Timeslab data sets to object names as described
#  in the Timeslab text.
#-----------------------------------------------------------------------
{
  air <- c(
   112, 118, 132, 129, 121, 135, 148, 148, 136, 119, 104, 118,
   115, 126, 141, 135, 125, 149, 170, 170, 158, 133, 114, 140,
   145, 150, 178, 163, 172, 178, 199, 199, 184, 162, 146, 166,
   171, 180, 193, 181, 183, 218, 230, 242, 209, 191, 172, 194,
   196, 196, 236, 235, 229, 243, 264, 272, 237, 211, 180, 201,
   204, 188, 235, 227, 234, 264, 302, 293, 259, 229, 203, 229,
   242, 233, 267, 269, 270, 315, 364, 347, 312, 274, 237, 278,
   284, 277, 317, 313, 318, 374, 413, 405, 355, 306, 271, 306,
   315, 301, 356, 348, 355, 422, 465, 467, 404, 347, 305, 336,
   340, 318, 362, 348, 363, 435, 491, 505, 404, 359, 310, 337,
   360, 342, 406, 396, 420, 472, 548, 559, 463, 407, 362, 405,
   417, 391, 419, 461, 472, 535, 622, 606, 508, 461, 390, 432)

artif <- c(
       12.5971,  12.2084,  13.0362,  10.1389,   8.8044,   7.5840,
        9.7160,  11.4201,  13.0817,  14.1539,  14.0661,  12.5436,
       10.7436,  10.2038,   9.7150,   9.7878,  10.0856,  11.2560,
       12.5911,  13.5058,  14.8224,  14.0889,  13.0398,  11.1967,
       10.3866,  10.0452,   9.1542,  12.0071,  13.9843,  14.1276,
       15.7731,  15.1213,  13.4067,  10.9186,   8.1324,   7.5325,
        8.2933,  11.5182,  13.5241,  14.6620,  17.4796,  17.6391,
       15.4382,  12.9454,  11.2046,  10.8988,  12.0149,  15.5140,
       19.5596,  20.1978,  18.5038,  19.2840,  16.5496,  15.9002,
       13.0752,  12.6914,  12.6856,  14.5315,  15.5767,  17.9975,
       19.2265,  18.5547,  19.3905,  16.4114,  15.6677,  15.4860,
       15.8087,  16.1202,  16.8718,  18.4903,  18.6455,  17.8001,
       16.9917,  14.7829,  15.3884,  15.4088,  15.1301,  15.2064,
       16.5537,  20.0116,  19.5788,  21.2045,  20.3659,  18.6175,
       17.0249,  18.1552,  18.8584,  18.8034,  18.7778,  19.9501,
       22.4343,  22.0911,  19.1892,  16.8242,  14.4648,  15.6066,
       15.4240,  17.3307,  20.0159,  21.5673)

bev <- c(
   17,  19,  20,  15,  13,  14,  14,  14,  14,  11,  16,  19,
   23,  18,  17,  20,  20,  18,  14,  16,  21,  24,  15,  16,
   20,  14,  16,25.5,25.8,  26,  26,  29,  20,  18,  16,  22,
   22,  16,  19,  17,  17,  19,  20,  24,  28,  36,  20,  14,
   18,  27,  29,  36,  29,  27,  30,  38,  50,  24,  25,  30,
   31,  37,  41,  36,  32,  47,  42,  37,  34,  36,  43,  55,
   64,  79,  59,  47,  48,  49,  45,  53,  55,  55,  54,  56,
   52,  76, 113,  68,  59,  74,  78,  69,  78,  73,  88,  98,
  109, 106,  87,  77,  77,  63,  70,  70,  63,  61,  66,  78,
   93,  97,  77,  83,  81,  82,  78,  75,  80,  87,  72,  65,
   74,  91, 115,  99,  99, 115, 101,  90,  95, 108, 147, 112,
  108,  99,  96, 102, 105, 114, 103,  98, 103, 101, 110, 109,
   98,  84,  90, 120, 124, 136, 120, 135, 100,  70,  60,  72,
   70,  71,  94,  95, 110, 154, 116,  99,  82,  76,  64,  63,
   68,  64,  67,  71,  72,  89, 114, 102,  85,  88,  97,  94,
   88,  79,  74,  79,  95,  70,  72,  63,  60,  74,  75,  91,
  126, 161, 109, 108, 110, 130, 166, 143, 103,  89,  76,  93,
   82,  71,  69,  75, 134, 183, 113, 108, 121, 139, 109,  90,
   88,  88,  93, 106,  89,  79,  91,  96, 111, 112, 104,  94,
   98,  88,  94,  81,  77,  84,  92,  96, 102,  95,  98, 125,
  162, 113,  94,  85,  89, 109, 110, 109, 120, 116, 101, 113,
  109, 105,  94, 102, 141, 135, 118, 115, 111, 127, 124, 113,
  122, 130, 137, 148, 142, 143, 176, 184, 164, 146, 147, 124,
  119, 135, 125, 116, 132, 133, 144, 145, 146, 138, 139, 154,
  181, 185, 151, 139, 157, 155, 191, 248, 185, 168, 176, 243,
  289, 251, 232, 207, 276, 250, 216, 205, 206, 208, 226, 302,
  261, 207, 209, 280, 381, 266, 197, 177, 170, 152, 156, 141,
  142, 137, 161, 189, 226, 194, 217, 199, 151, 144, 138, 145,
  156, 184, 216, 204, 186, 197, 183, 175, 183, 230, 278, 179,
  161, 150, 159, 180, 223, 294, 300, 297, 232, 179, 180, 215,
  258, 236, 202, 174, 179, 210, 268, 267, 208, 224)

buffsnow <- c(
    126.4,  82.4,  78.1,  51.1,  90.9,  76.2, 104.5,  87.4, 110.5,  25.0,
     69.3,  53.5,  39.8,  63.6,  46.7,  72.9,  79.6,  83.6,  80.7,  60.3,
     79.0,  74.4,  49.6,  54.7,  71.8,  49.1, 103.9,  51.6,  82.4,  83.6,
     77.8,  79.3,  89.6,  85.5,  58.0, 120.7, 110.5,  65.4,  39.9,  40.1,
     88.7,  71.4,  83.0,  55.9,  89.9,  84.8, 105.2, 113.7, 124.7, 114.5,
    115.6, 102.4, 101.4,  89.8,  71.5,  70.9,  98.3,  55.5,  66.1,  78.4,
    120.5,  97.0, 110.0)

calfem <- c(
 35., 32., 30., 31., 44., 29., 45., 43., 38., 27.,
 38., 33., 55., 47., 45., 37., 50., 43., 41., 52.,
 34., 53., 39., 32., 37., 43., 39., 35., 44., 38.,
 24., 23., 31., 44., 38., 50., 38., 51., 31., 31.,
 51., 36., 45., 51., 34., 52., 47., 45., 46., 39.,
 48., 37., 35., 52., 42., 45., 39., 37., 30., 35.,
 28., 45., 34., 36., 50., 44., 39., 32., 39., 45.,
 43., 39., 31., 27., 30., 42., 46., 41., 36., 45.,
 46., 43., 38., 34., 35., 56., 36., 32., 50., 41.,
 39., 41., 47., 34., 36., 33., 35., 38., 38., 34.,
 53., 34., 34., 38., 35., 32., 42., 34., 46., 30.,
 46., 45., 54., 34., 37., 35., 40., 42., 58., 51.,
 32., 35., 38., 33., 39., 47., 38., 52., 30., 34.,
 40., 35., 42., 41., 42., 38., 24., 34., 43., 36.,
 55., 41., 45., 41., 37., 43., 39., 33., 43., 40.,
 38., 45., 46., 34., 35., 48., 51., 36., 33., 46.,
 42., 48., 34., 41., 35., 40., 34., 30., 36., 40.,
 39., 45., 38., 47., 33., 30., 42., 43., 41., 41.,
 59., 43., 45., 38., 37., 45., 42., 57., 46., 51.,
 41., 47., 26., 35., 44., 41., 42., 36., 45., 45.,
 45., 47., 38., 42., 35., 36., 39., 45., 43., 47.,
 36., 41., 50., 39., 41., 46., 64., 45., 34., 38.,
 44., 48., 46., 44., 37., 39., 44., 45., 33., 44.,
 38., 46., 46., 40., 39., 44., 48., 50., 41., 42.,
 51., 41., 44., 38., 68., 40., 42., 51., 44., 45.,
 36., 57., 44., 42., 53., 42., 34., 40., 56., 44.,
 53., 55., 39., 59., 55., 73., 55., 44., 43., 40.,
 47., 51., 56., 49., 54., 56., 47., 44., 43., 42.,
 45., 50., 48., 43., 40., 59., 41., 42., 51., 49.,
 45., 43., 42., 38., 47., 38., 36., 42., 35., 28.,
 44., 36., 45., 46., 48., 49., 43., 42., 59., 45.,
 52., 46., 42., 40., 40., 45., 35., 35., 40., 39.,
 33., 42., 47., 51., 44., 40., 57., 49., 45., 49.,
 51., 46., 44., 52., 45., 32., 46., 41., 34., 33.,
 36., 49., 43., 43., 34., 39., 35., 52., 47., 52.,
 39., 40., 42., 42., 53., 39., 40., 38., 44., 34.,
 37., 52., 48., 55., 50)

cos4 <- c(
     16.0000,   8.5679,  -3.5343,   6.3660,    .3906,  -2.1213,
     14.0981,   3.7592,  -9.5681,   1.0000,  -2.7192,  -2.1210,
     17.0981,   8.7175,  -4.3377,   4.6340,  -2.2175,  -5.5260,
     10.0000,   -.9077, -14.6619,  -4.3660,  -8.1944,  -7.5390,
     11.9019,   3.9011,  -8.6281,   1.0000,  -5.0847,  -7.5392,
      8.9019,  -1.0573, -13.8585,  -2.6339,  -5.5864,  -4.1342,
     16.0000,   8.5679,  -3.5343,   6.3660,    .3906,  -2.1212,
     14.0981,   3.7592,  -9.5681,   1.0000,  -2.7192,  -2.1210,
     17.0981,   8.7175,  -4.3377,   4.6340,  -2.2175,  -5.5260,
     10.0000,   -.9077, -14.6619,  -4.3660,  -8.1944,  -7.5390,
     11.9020,   3.9011,  -8.6281,   1.0000,  -5.0847,  -7.5392,
      8.9020,  -1.0573, -13.8585,  -2.6340,  -5.5864,  -4.1342,
     16.0000,   8.5678,  -3.5342,   6.3660,    .3906,  -2.1212,
     14.0981,   3.7592,  -9.5681,   1.0000,  -2.7192,  -2.1210,
     17.0981,   8.7175,  -4.3377,   4.6340,  -2.2175,  -5.5260,
     10.0000,   -.9078, -14.6619,  -4.3660,  -8.1944,  -7.5389,
     11.9020,   3.9011,  -8.6281,   1.0000,  -5.0847,  -7.5392,
      8.9020,  -1.0573, -13.8585,  -2.6339,  -5.5864,  -4.1342,
     16.0000,   8.5678,  -3.5342,   6.3660,    .3905,  -2.1212,
     14.0981,   3.7592,  -9.5681,   1.0000,  -2.7192,  -2.1210,
     17.0981,   8.7175,  -4.3377,   4.6340,  -2.2174,  -5.5259,
     10.0001,   -.9078, -14.6619,  -4.3660,  -8.1944,  -7.5390,
     11.9020,   3.9009,  -8.6281,   1.0000,  -5.0847,  -7.5392,
      8.9019,  -1.0575, -13.8584,  -2.6339,  -5.5864,  -4.1342)



cradfq <- c(
      5.40,  5.40,  5.10,  4.80,  5.50,  6.10,  6.40,  6.80,
      6.40,  6.00,  6.10,  6.10,  5.30,  5.50,  5.70,  5.00,
      5.80,  9.00, 10.40,  9.90, 10.70,  9.70,  8.80,  8.10,
      8.40,  6.40,  6.70,  7.40,  9.10, 11.90, 12.50, 11.80,
     11.10, 13.10, 12.30,  9.70,  7.70,  7.60,  8.10,  8.90,
     10.40, 13.00, 13.20, 12.30, 11.50, 12.80, 12.20, 10.00,
      8.10,  7.30,  7.40,  7.90,  9.60, 12.10, 13.10, 12.40,
     11.10, 11.70, 10.30,  9.30,  8.00,  7.40,  6.90,  7.60,
      8.90, 12.20, 11.40, 10.30,  9.50,  9.60,  9.80,  8.40,
      7.40,  6.70,  6.40,  7.10,  9.20, 10.70, 10.60,  9.30,
      9.00,  8.90,  8.10,  7.20,  5.70,  5.80,  5.80,  6.20,
      6.90,  9.10,  9.00,  8.70,  8.00,  7.80,  8.70,  7.40,
      6.10,  5.30,  4.90,  5.10,  6.00,  7.40,  7.50,  7.00,
      6.40,  6.60,  6.20,  5.60,  5.10,  5.20,  5.20,  4.80,
      5.20,  6.00,  6.60,  7.10,  6.40,  6.00,  5.90,  5.20,
      4.90,  4.70,  4.70,  5.00,  5.60,  7.00,  6.60,  7.10,
      6.50,  6.60,  6.20,  6.30,  5.70,  5.80,  5.50,  5.70,
      6.40,  9.30,  9.30,  8.60,  7.70,  9.40,  9.70,  9.20,
      6.20,  6.30,  6.10,  7.20,  8.60, 11.50, 12.20, 12.00,
     12.00, 13.00, 12.50, 11.20, 10.00,  6.90,  7.10,  7.70,
     10.60, 13.00, 13.50, 12.60, 12.30, 11.30, 11.20, 10.80,
      8.20,  7.60,  7.10,  6.80,  9.20, 11.50, 12.40, 11.70,
     12.10, 12.50, 12.50, 10.60,  8.30,  7.00,  7.10,  7.80,
     10.10, 12.60, 13.10, 12.00, 11.30, 11.00, 10.50,  8.80,
      7.20,  6.80,  6.00,  6.10,  7.00,  9.20,  9.00,  8.60,
      8.40,  8.70,  7.80,  6.30,  6.10,  6.00,  5.60,  5.80,
      7.00,  8.60,  9.20,  8.80,  8.00,  7.50,  6.50,  5.50,
      4.80,  5.20,  5.00,  5.30,  6.00,  7.20,  7.60,  7.30,
      6.40,  6.20,  5.20,  5.20,  4.80,  4.80,  4.30,  4.70,
      6.00,  6.20,  6.50,  6.30,  5.90,  5.70,  5.50,  5.00)

cradfqsun <- c(
      15.9850,  15.3350,  15.4650,  15.4150,  15.2000,  15.2850,
      15.4350,  15.7700,  15.9300,  16.0250,  16.1550,  15.6100,
      16.3650,  17.2850,  16.6950,  16.5050,  17.1050,  17.6600,
      18.2100,  18.0750,  18.1400,  18.7150,  18.8550,  18.7450,
      17.7300,  18.5000,  17.6150,  19.3500,  18.8000,  19.4500,
      20.7700,  21.1700,  21.6250,  21.4250,  19.1950,  20.4650,
      20.8350,  21.5150,  22.2550,  21.8850,  20.0350,  21.2450,
      18.7200,  19.4400,  19.9200,  20.9600,  19.3250,  20.0500,
      21.3700,  19.8750,  23.2650,  20.7850,  19.4800,  19.9550,
      21.1100,  19.6350,  19.0150,  18.8700,  18.2300,  20.4550,
      20.9150,  20.0500,  19.8800,  20.2900,  20.6300,  19.4050,
      18.4050,  17.1050,  17.5250,  17.9700,  19.1650,  18.0350,
      17.7200,  19.1950,  18.3750,  20.2750,  18.3250,  17.7500,
      17.9200,  18.4150,  17.2800,  17.2250,  17.3200,  16.6400,
      16.4750,  17.9900,  18.3450,  18.0000,  18.2950,  17.3150,
      16.9200,  16.6850,  16.7800,  17.6400,  17.7100,  18.0350,
      16.2500,  15.5700,  15.8850,  16.0100,  15.8600,  15.9600,
      16.5350,  16.1250,  15.6200,  16.4450,  16.3700,  16.3050,
      15.7050,  15.3800,  15.6600,  15.9700,  15.5000,  15.3900,
      15.5100,  15.9400,  15.1850,  15.0250,  15.5500,  15.0150,
      15.1250,  15.2500,  15.2500,  15.8350,  15.7150,  15.8450,
      15.5400,  16.4200,  15.9250,  15.6350,  16.0750,  16.6000,
      16.5300,  16.8100,  17.1300,  16.2950,  16.7450,  18.4400,
      17.3000,  16.3700,  17.3800,  19.3100,  18.8300,  18.7850,
      19.2450,  18.6750,  20.8100,  20.3600,  19.7200,  20.1150,
      21.1900,  21.0850,  20.7850,  21.6700,  21.4900,  22.4900,
      25.0650,  23.1950,  22.8950,  24.4400,  23.4700,  23.1800,
      21.4000,  20.8250,  20.4250,  19.3050,  19.7400,  24.4850,
      23.7000,  23.3900,  22.1100,  22.8950,  22.1650,  21.8150,
      19.7900,  21.9000,  20.9550,  24.1150,  22.8750,  22.3500,
      20.3100,  21.0850,  21.2900,  21.1900,  22.2650,  21.5800,
      22.1750,  20.8800,  20.0800,  19.7400,  20.4850,  20.6700,
      20.3100,  19.1800,  19.5500,  19.2600,  17.5650,  18.0700,
      17.7400,  17.7050,  17.9950,  17.9950,  17.7950,  19.6450,
      20.4250,  20.0300,  18.0750,  18.0500,  19.1550,  17.5800,
      17.6200,  17.2900,  17.0350,  16.1350,  16.1000,  16.4550,
      16.1700,  16.8200,  16.9650,  17.7450,  16.4100,  16.1900,
      16.1050,  16.7150,  16.3250,  15.1950,  15.5000,  16.3900,
      15.6250,  16.0900,  15.4300,  16.1750,  15.9650,  15.4100,
      15.0800,  15.1250,  15.0100,  15.0250,  15.5450,  15.0900)


eriel <- c(
    14.763,14.649,15.085,16.376,16.926,16.774,16.490,15.769,15.180,14.383,
    14.478,14.364,13.928,13.283,13.700,15.465,16.243,16.490,16.243,15.787,
    15.446,14.649,13.776,13.188,13.283,12.657,12.979,13.909,14.535,14.877,
    14.858,14.288,13.947,13.416,12.903,13.454,13.491,13.567,13.397,14.440,
    15.161,15.427,15.693,15.142,14.763,14.288,13.074,12.543,12.239,12.011,
    12.827,13.567,13.548,13.302,13.188,13.112,12.827,12.201,11.917,11.803,
    11.157,10.892,11.120,12.600,13.283,13.416,13.340,13.529,13.776,14.307,
    13.852,13.833,13.169,12.941,13.188,14.383,14.763,15.218,15.161,14.858,
    14.156,13.586,13.150,14.137,14.231,14.364,13.833,14.478,15.009,15.617,
    16.148,15.977,15.142,14.592,14.364,14.497,14.554,14.991,15.863,17.932,
    19.184,19.184,18.956,18.254,17.514,16.660,16.338,16.319,17.457,17.173,
    17.856,18.596,18.558,18.159,17.685,16.812,16.072,15.332,14.478,14.213,
    13.738,13.169,12.581,13.245,13.852,14.175,14.288,13.985,13.435,12.884,
    12.429,12.410,13.397,13.909,13.833,14.099,14.687,14.611,14.383,13.909,
    13.359,12.296,12.106,11.803,12.353,12.220,12.827,14.250,15.085,14.953,
    14.440,13.890,13.036,12.201,11.404,11.309,10.987,10.361,10.304,11.347,
    11.784,11.841,11.841,11.651,11.404,10.873,10.209,10.076,10.247,10.133,
    10.740,11.556,12.201,12.505,12.732,12.201,12.068,11.290,11.139,11.101,
    10.342,10.000,11.347,12.770,13.321,13.340,13.188,12.676,12.315,12.049,
    11.594,11.252,12.467,13.491,13.491,14.156,15.256,15.598,16.034,15.598,
    14.516,13.435,12.638,12.106,12.144,12.979,11.917,15.066,15.199,15.427,
    15.427,15.408,14.706,14.137,13.302,12.979,13.036,12.903,13.719,14.801,
    15.541,15.617,15.503,15.218,14.592,13.852,13.416,13.055,12.315,12.239,
    12.562,13.890,14.687,15.313,15.408,15.028,14.782,14.213,13.435,13.662,
    14.288,13.491,13.150,13.586,13.871,14.175,14.099,13.719,13.093,12.562,
    12.030,12.144,11.860,12.410,12.770,14.630,15.218,15.920,15.882,15.750,
    15.275,14.801,14.554,14.345,15.104,14.554,14.953,15.958,17.590,18.805,
    18.805,18.330,17.476,16.812,16.224,15.541,14.744,14.478,14.725,16.319,
    17.324,17.609,17.211,16.546,15.977,15.427,14.972,14.535,14.213,13.719,
    15.009,16.319,17.078,17.913,18.159,17.704,17.059,17.268,16.546,16.072,
    15.996,15.294,15.901,16.300,16.546,17.495,17.666,16.964,16.148,15.427,
    14.839,14.269,14.118,14.156,14.080,16.414,18.216,19.298,18.824,18.368,
    17.875,16.869,16.129,15.712,15.617,15.161,16.148,17.609,18.406,18.463,
    18.254,17.609,16.812,15.712,15.066,14.763,14.877,15.370,15.731,15.996,
    16.224,16.186,16.015,15.446,14.573,14.080,13.226,12.979,14.459,15.655,
    15.636,17.154,17.381,17.116,16.736,16.167,15.920,15.408,15.066,15.769,
    15.787,15.863,17.154,18.178,18.653,18.615,18.406,17.837,17.078,16.471,
    16.224,16.395,17.400,18.672,19.089,19.829,20.000,19.943,19.526,18.975,
    18.463,17.268,16.414,16.414,16.755,16.736,17.173,17.647,18.216,18.767,
    18.539,18.273,17.419,16.679,15.996,15.465,15.408,15.123,16.072,17.685,
    18.235,18.064,17.818,17.438,16.812,17.116,17.211,17.097,17.495,17.097,
    18.121,18.748,18.767,18.539,18.102,17.818,17.002,16.357,15.560,15.313,
    14.782,13.719,14.839,15.825,17.324,17.723,17.685,17.514,16.907,15.863,
    14.934,14.706,14.383,14.345,14.763,15.958,16.698,16.793,17.230,16.509,
    15.787,14.991,14.080,14.326,14.497,13.510,13.662,14.042,14.269,14.478,
    14.972,14.972,14.573,13.776,13.017,12.600,12.296,13.131,13.966,15.028,
    15.844,15.769,15.237,14.801,14.137,13.890,13.416,13.871,14.478,14.725,
    14.763,15.806,16.565,17.097,17.306,17.211,16.717,15.787,14.801,14.231,
    13.966,14.004,15.313,16.357,17.742,17.609,17.249,17.078,16.509,15.465,
    14.706,14.213,13.662,13.928,14.516,15.180,15.351,15.579,15.446,15.199,
    14.725,14.383,14.231,13.776,13.150,12.713,13.283,14.478,14.858,14.782,
    14.307,14.023,13.529,12.922,12.410,11.936,11.784,11.992,12.619,13.662,
    14.231,14.099,13.833,13.416,12.922,11.974,11.461,11.423,11.860,12.163,
    13.226,13.966,14.535,14.516,14.231,13.966,13.738,13.226,12.998,13.131,
    13.700,13.814,14.383,15.047,15.693,15.844,15.712,15.313,14.763,13.548,
    13.586,14.516,14.383,14.497,14.744,15.806,16.527,16.546,16.717,16.433,
    15.769,15.275,15.123,15.541,15.825,16.376,16.357,16.907,17.021,17.419,
    17.533,17.268,16.603,15.825,15.446,15.636,15.693,16.755,16.509,17.723,
    18.691,19.127,19.564,19.203,18.216,17.211,16.660,16.831,15.769,15.731,
    15.996,17.021,17.552,17.837,17.856,17.571,17.078,16.660,16.433,16.584)

gasfurnin <- c(
     -.11,    .00,    .18,    .34,    .37,    .44,    .46,    .35,
      .13,   -.18,   -.59,  -1.05,  -1.42,  -1.52,  -1.30,   -.81,
     -.47,   -.19,    .09,    .44,    .77,    .87,    .88,    .89,
      .99,   1.26,   1.77,   1.98,   1.93,   1.87,   1.83,   1.77,
     1.61,   1.26,    .79,    .36,    .12,    .09,    .33,    .64,
      .96,   1.41,   2.67,   2.83,   2.81,   2.48,   1.93,   1.49,
     1.21,   1.24,   1.61,   1.90,   2.02,   1.82,    .54,    .12,
      .01,    .16,    .67,   1.02,   1.15,   1.15,   1.11,   1.12,
     1.22,   1.26,   1.16,    .91,    .62,    .25,   -.28,  -1.08,
    -1.55,  -1.80,  -1.83,  -1.46,   -.94,   -.57,   -.43,   -.58,
     -.96,  -1.62,  -1.88,  -1.89,  -1.75,  -1.47,  -1.20,   -.93,
     -.52,    .04,    .79,    .94,    .93,   1.01,   1.14,   1.20,
     1.05,    .60,   -.08,   -.31,   -.29,   -.15,   -.11,   -.19,
     -.25,   -.23,   -.01,    .25,    .33,    .10,   -.42,  -1.14,
    -2.28,  -2.59,  -2.72,  -2.51,  -1.79,  -1.35,  -1.08,   -.91,
     -.88,   -.88,   -.80,   -.54,   -.42,   -.27,    .00,    .40,
      .84,   1.28,   1.61,   1.75,   1.68,   1.49,    .99,    .65,
      .58,    .58,    .63,    .75,    .90,    .99,    .97,    .79,
      .40,   -.16,   -.55,   -.60,   -.42,   -.19,   -.05,    .06,
      .16,    .30,    .52,    .57,    .56,    .57,    .59,    .67,
      .93,   1.34,   1.46,   1.35,    .77,    .22,   -.24,   -.71,
    -1.10,  -1.27,  -1.17,   -.68,    .03,    .56,    .64,    .48,
      .11,   -.31,   -.70,  -1.05,  -1.22,  -1.18,   -.87,   -.34,
      .06,    .08,    .00,    .00,    .21,    .56,    .78,    .86,
      .92,    .86,    .42,   -.34,   -.96,  -1.81,  -2.38,  -2.50,
    -2.47,  -2.33,  -2.05,  -1.74,  -1.26,   -.57,   -.14,   -.02,
     -.05,   -.14,   -.28,   -.53,   -.87,  -1.24,  -1.44,  -1.42,
    -1.17,   -.81,   -.63,   -.58,   -.63,   -.71,   -.85,  -1.04,
    -1.35,  -1.63,  -1.62,  -1.15,   -.49,   -.16,   -.01,   -.09,
     -.62,  -1.09,  -1.52,  -1.86,  -2.03,  -2.02,  -1.96,  -1.95,
    -1.79,  -1.30,  -1.03,   -.92,   -.80,   -.87,  -1.05,  -1.12,
     -.88,   -.40,    .19,    .66,    .71,    .61,    .50,    .60,
      .94,   1.22,   1.25,    .82,    .10,    .03,    .38,    .92,
     1.03,    .87,    .53,    .09,   -.46,   -.75,   -.95,  -1.03,
     -.93,   -.64,   -.42,   -.28,   -.16,   -.03,    .10,    .25,
      .28,    .00,   -.49,   -.76,   -.82,   -.74,   -.53,   -.20,
      .03,    .20,    .25,    .19,    .13,    .02,   -.18,   -.26)

gasfurnout <- c(
53.8, 53.6, 53.5, 53.5, 53.4, 53.1, 52.7, 52.4, 52.2, 52.0, 52.0, 52.4,
53.0, 54.0, 54.9, 56.0, 56.8, 56.8, 56.4, 55.7, 55.0, 54.3, 53.2, 52.3,
51.6, 51.2, 50.8, 50.5, 50.0, 49.2, 48.4, 47.9, 47.6, 47.5, 47.5, 47.6,
48.1, 49.0, 50.0, 51.1, 51.8, 51.9, 51.7, 51.2, 50.0, 48.3, 47.0, 45.8,
45.6, 46.0, 46.9, 47.8, 48.2, 48.3, 47.9, 47.2, 47.2, 48.1, 49.4, 50.6,
51.5, 51.6, 51.2, 50.5, 50.1, 49.8, 49.6, 49.4, 49.3, 49.2, 49.3, 49.7,
50.3, 51.3, 52.8, 54.4, 56.0, 56.9, 57.5, 57.3, 56.6, 56.0, 55.4, 55.4,
56.4, 57.2, 58.0, 58.4, 58.4, 58.1, 57.7, 57.0, 56.0, 54.7, 53.2, 52.1,
51.6, 51.0, 50.5, 50.4, 51.0, 51.8, 52.4, 53.0, 53.4, 53.6, 53.7, 53.8,
53.8, 53.8, 53.3, 53.0, 52.9, 53.4, 54.6, 56.4, 58.0, 59.4, 60.2, 60.0,
59.4, 58.4, 57.6, 56.9, 56.4, 56.0, 55.7, 55.3, 55.0, 54.4, 53.7, 52.8,
51.6, 50.6, 49.4, 48.8, 48.5, 48.7, 49.2, 49.8, 50.4, 50.7, 50.9, 50.7,
50.5, 50.4, 50.2, 50.4, 51.2, 52.3, 53.2, 53.9, 54.1, 54.0, 53.6, 53.2,
53.0, 52.8, 52.3, 51.9, 51.6, 51.6, 51.4, 51.2, 50.7, 50.0, 49.4, 49.3,
49.7, 50.6, 51.8, 53.0, 54.0, 55.3, 55.9, 55.9, 54.6, 53.5, 52.4, 52.1,
52.3, 53.0, 53.8, 54.6, 55.4, 55.9, 55.9, 55.2, 54.4, 53.7, 53.6, 53.6,
53.2, 52.5, 52.0, 51.4, 51.0, 50.9, 52.4, 53.5, 55.6, 58.0, 59.5, 60.0,
60.4, 60.5, 60.2, 59.7, 59.0, 57.6, 56.4, 55.2, 54.5, 54.1, 54.1, 54.4,
55.5, 56.2, 57.0, 57.3, 57.4, 57.0, 56.4, 55.9, 55.5, 55.3, 55.2, 55.4,
56.0, 56.5, 57.1, 57.3, 56.8, 55.6, 55.0, 54.1, 54.3, 55.3, 56.4, 57.2,
57.8, 58.3, 58.6, 58.8, 58.8, 58.6, 58.0, 57.4, 57.0, 56.4, 56.3, 56.4,
56.4, 56.0, 55.2, 54.0, 53.0, 52.0, 51.6, 51.6, 51.1, 50.4, 50.0, 50.0,
52.0, 54.0, 55.1, 54.5, 52.8, 51.4, 50.8, 51.2, 52.0, 52.8, 53.8, 54.5,
54.9, 54.9, 54.8, 54.4, 53.7, 53.3, 52.8, 52.6, 52.6, 53.0, 54.3, 56.0,
57.0, 58.0, 58.6, 58.5, 58.3, 57.8, 57.3, 57.0)

lh1 <- c(
 14, 15, 16, 14, 12, 10, 11, 16, 13, 12, 19, 19, 14, 13, 14, 20, 16,
 13, 13, 20, 16, 14, 13, 19, 18, 12, 13, 20, 17, 18, 13, 20, 14, 13,
 11, 11, 18, 14, 11,  9, 20, 15, 14, 11, 18, 13, 13,  9, 19, 14, 12,
 10, 17, 16, 12, 11, 10, 19, 16, 12, 11, 17, 15, 13, 13, 19, 17, 11,
 11, 17, 14, 11,  9,  9, 14, 11,  9,  8, 12, 11,  9,  9, 13, 11,  8,
  9, 16, 13, 11, 10, 13, 13, 10,  9, 12, 15, 11, 11, 13, 13, 11, 11,
  9, 12, 12, 10, 10, 16, 13, 12, 10,  9, 12, 11, 10,  9, 10, 14, 13,
 11, 10,  8, 16, 12, 12, 10,  9,  9, 12, 12, 10,  9, 10,  9, 12, 15,
 12, 11,  9, 10,  8,  8, 12, 12)

lh2 <- c(
 14, 11, 11, 16, 13, 12, 10,  9, 20, 16, 13, 13, 11,  9, 16, 15, 13,
 14, 10, 12, 18, 13, 13, 11, 12, 20, 20, 16, 14, 13, 21, 17, 14, 14,
 12, 19, 18, 16, 15, 14, 16, 20, 17, 14, 15, 12, 18, 18, 15, 13, 13,
 22, 16, 15, 14, 12, 20, 19, 16, 14, 13, 12, 21, 19, 18, 15, 15, 25,
 21, 17, 14, 18, 19, 14, 14, 14, 19, 19, 17, 14, 13, 21, 17, 16, 14,
 16, 19, 18, 14, 13, 18, 18, 15, 16, 12, 19, 16, 14, 12, 18, 18, 15,
 16, 12, 14, 17, 15, 13, 11, 20, 17, 14, 13, 14, 14, 14, 16, 13, 16,
 17, 15, 13, 11, 16, 19, 16, 13, 14, 18, 17, 16, 15, 14, 19, 19, 15,
 13, 13, 11, 17, 17, 15, 13, 19)

lynx <- c(
      269.,     321.,     585.,     871.,    1475.,    2821.,
     3928.,    5943.,    4950.,    2577.,     523.,      98.,
      184.,     279.,     409.,    2285.,    2685.,    3409.,
     1824.,     409.,     151.,      45.,      68.,     213.,
      546.,    1033.,    2129.,    2536.,     957.,     361.,
      377.,     225.,     360.,     731.,    1638.,    2725.,
     2871.,    2119.,     684.,     299.,     236.,     245.,
      552.,    1623.,    3311.,    6721.,    4245.,     687.,
      255.,     473.,     358.,     784.,    1594.,    1676.,
     2251.,    1426.,     756.,     299.,     201.,     229.,
      469.,     736.,    2042.,    2811.,    4431.,    2511.,
      389.,      73.,      39.,      49.,      59.,     188.,
      377.,    1292.,    4031.,    3495.,     587.,     105.,
      153.,     387.,     758.,    1307.,    3465.,    6991.,
     6313.,    3794.,    1836.,     345.,     382.,     808.,
     1388.,    2713.,    3800.,    3091.,    2985.,    3790.,
      674.,      81.,      80.,     108.,     229.,     399.,
     1132.,    2432.,    3574.,    2935.,    1537.,     529.,
      485.,     662.,    1000.,    1590.,    2657.,    3396)

mlco2 <- c(
   14.88, 15.62, 16.33, 17.59, 17.93, 17.71, 15.92, 15.15, 14.02, 12.83,
   13.64, 14.71, 15.62, 16.59, 16.94, 17.77, 18.29, 18.24, 16.67, 14.96,
   14.12, 13.58, 15.14, 15.77, 16.62, 17.16, 17.90, 19.21, 20.02, 19.74,
   18.15, 16.00, 14.23, 14.07, 15.04, 16.19, 16.97, 17.74, 18.63, 19.43,
   20.47, 19.71, 18.78, 16.84, 15.16, 15.56, 16.14, 17.13, 18.06, 18.59,
   19.74, 20.63, 21.21, 20.83, 19.55, 17.75, 16.27, 15.62, 16.84, 17.70,
   18.80, 19.08, 20.15, 21.49, 22.25, 21.50, 19.67, 17.61, 16.25, 16.17,
   17.01, 18.36, 19.37, 19.93, 20.40, 21.65, 22.26, 22.19, 20.49, 18.48,
   17.13, 17.02, 17.84, 18.78, 19.55, 20.65, 21.15, 22.31, 22.35, 22.19,
   21.53, 19.13, 17.99, 17.70, 19.15, 19.27, 20.22, 21.23, 22.13, 23.30,
   23.57, 23.29, 22.36, 19.71, 17.89, 17.54, 19.36, 20.51, 21.80, 22.03,
   22.50, 24.00, 24.46, 23.46, 22.19, 20.57, 18.91, 18.81, 20.24, 21.59,
   22.15, 22.73, 23.50, 24.52, 25.11, 25.06, 23.62, 21.55, 19.89, 19.80,
   20.73, 22.25, 23.73, 24.53, 25.62, 26.58, 27.24, 26.53, 25.63, 23.28,
   22.21, 21.67, 22.61, 24.07, 24.91, 25.81, 26.85, 28.07, 27.97, 27.77,
   26.44, 24.92, 23.49, 23.50, 24.34, 25.39, 26.46, 26.93, 27.56, 28.23,
   29.51, 29.04, 27.87, 26.00, 24.06, 24.20, 25.48, 26.62, 27.30, 28.20,
   28.50, 30.22, 30.58, 29.48, 28.56, 26.77, 25.39, 25.72, 26.97, 28.09,
   29.16, 30.02, 30.95, 31.95, 32.85, 32.58, 31.30, 29.64, 28.12, 27.67,
   28.69, 29.07, 29.84, 31.17, 31.94, 33.20, 33.55, 32.76, 31.79, 29.61,
   27.88, 27.85, 28.79, 30.13, 30.64, 31.21, 31.88, 33.13, 33.76, 33.74,
   32.05, 30.40, 28.89, 28.55, 29.57, 30.83)

normwn <- c(
        .0943,   -.2142,    .8439,    .7062,    .5208,    .2100,
       1.4550,    .4649,   -.4410,   -.9391,   -.5790,  -1.0580,
       -.6708,    .8435,  -1.7748,   -.8621,   -.0782,   -.2264,
       -.5232,   -.7443,    .4134,   -.3300,    .2164,   -.1358,
       -.6790,    .8799,  -1.4526,    .9732,    .0566,  -1.4173,
        .6477,   -.2606,   -.4727,   -.0500,    .8909,    .6899,
      -1.0235,   1.2896,   -.2868,   -.7879,   1.7034,    .7273,
       -.5479,   -.6860,    .7058,    .2904,    .2679,   1.8345,
       2.5018,   -.0287,  -1.6743,   1.3541,   -.9202,   1.1854,
      -1.0682,    .0897,   -.6558,    .0990,   -.9942,    .5082,
        .4487,   -.1672,   2.2177,   -.6841,    .3455,   -.2513,
        .1213,   -.9045,  -1.1407,   -.2203,   -.7520,   -.7361,
        .1348,   -.7635,   1.5830,    .5722,   -.7494,  -1.5595,
       1.0255,   -.0823,  -1.2096,   1.7183,  -1.0408,  -1.0837,
        .0935,   1.8412,    .6691,  -1.0086,  -1.7471,   -.6900,
       1.4420,    .2133,  -1.3858,   -.8704,  -1.3716,    .9627,
      -1.1039,   -.1556,   1.1266,    .4647,   -.3326,   -.0414,
       2.1232,   -.6909,    .2319,   -.5646)

nyct <- c(
    11.506,11.022,14.405,14.442,16.524,17.918,18.959,18.309,
    18.160,16.691,14.480,17.862,12.082,10.558,12.138,14.442,
    16.152,17.714,19.015,19.164,17.900,16.933,13.364,11.468,
    10.000,10.985,12.993,14.480,16.134,17.862,19.238,19.089,
    18.067,15.576,14.926,12.398,12.435,12.416,13.067,15.093,
    16.766,18.680,19.796,19.424,17.398,16.877,13.866,12.658,
    13.048,11.301,12.063,14.164,16.041,18.067,19.108,18.736,
    17.212,16.394,13.922,11.691,11.970,11.970,12.900,14.888,
    16.747,17.955,19.201,18.903,17.714,16.097,13.234,12.435,
    11.952,11.970,12.639,15.204,16.283,18.699,19.851,18.959,
    18.030,15.316,14.126,12.323,12.230,12.230,13.234,14.647,
    16.729,18.494,19.312,19.164,18.011,16.338,14.368,12.751,
    11.022,12.658,12.844,15.037,16.190,18.253,19.182,18.513,
    17.621,16.506,13.717,11.840,10.855,11.636,12.825,15.019,
    16.952,17.937,20.000,19.535,17.770,16.283,13.327,10.818,
    11.264,12.026,12.230,14.126,15.818,18.420,18.699,18.903,
    17.175,16.041,13.866,12.844,10.613,12.082,12.955,15.000,
    16.747,18.792,19.312,18.680,18.067,15.632,14.349,12.732,
    11.152,10.372,12.639,14.981,16.097,17.602,19.331,18.978,
    17.695,15.539,14.182,10.781,11.059,11.152,12.621,14.981,
    17.100,18.160,19.015,19.331,18.346,16.115,13.532,12.398)

nycb <- c(
    26.663,23.598,26.931,24.740,25.806,24.364,24.477,23.901,
    23.175,23.227,21.672,21.870,21.439,21.089,23.709,21.669,
    21.752,20.761,23.479,23.824,23.105,23.110,21.759,22.073,
    21.937,20.035,23.590,21.672,22.222,22.123,23.950,23.504,
    22.238,23.142,21.059,21.573,21.548,20.000,22.424,20.615,
    21.761,22.874,24.104,23.748,23.262,22.907,21.519,22.025,
    22.604,20.894,24.677,23.673,25.320,23.583,24.671,24.454,
    24.122,24.252,22.084,22.991,23.287,23.049,25.076,24.037,
    24.430,24.667,26.451,25.618,25.014,25.110,22.964,23.981,
    23.798,22.270,24.775,22.646,23.988,24.737,26.276,25.816,
    25.210,25.199,23.162,24.707,24.364,22.644,25.565,24.062,
    25.431,24.635,27.009,26.606,26.268,26.462,25.246,25.180,
    24.657,23.304,26.982,26.199,27.210,26.122,26.706,26.878,
    26.152,26.379,24.712,25.688,24.990,24.239,26.721,23.475,
    24.767,26.219,28.361,28.599,27.914,27.784,25.693,26.881,
    26.217,24.218,27.914,26.975,28.527,27.139,28.982,28.169,
    28.056,29.136,26.291,26.987,26.589,24.848,27.543,26.896,
    28.878,27.390,28.065,28.141,29.048,28.484,26.634,27.735,
    27.132,24.924,28.963,26.589,27.931,28.009,29.229,28.759,
    28.405,27.945,25.912,26.619,26.076,25.286,27.660,25.951,
    26.398,25.565,28.865,30.000,29.261,29.012,26.992,27.897)

raineast <- c(
  22.54,-17.46,-26.46, 11.54,  6.54,-26.46,  9.54, -1.46,-33.46,-14.46,
   6.54, -9.46,  6.54, -1.46, -1.46,   .54, -6.46,-15.46, -5.46,  4.54,
  -7.46, -2.46, -6.46, -2.46,  6.54,  1.54,  7.54, -9.46, -6.46,  9.54,
  10.54, -2.46,  -.46, 14.54,-13.46,  1.54, -4.46, -5.46,   .54,-13.46,
   5.54,  6.54, 14.54, -7.46, -3.46,  3.54,  3.54,-14.46, 16.54,  -.46,
   2.54, 12.54,  6.54, 13.54,  1.54,  -.46, 16.54,  5.54,  7.54, 15.54,
   7.54, 15.54, -2.46,  4.54,  3.54,  6.54,  9.54,  6.54,  3.54,  2.54,
  -2.46,  7.54,  8.54, 10.54,   .54,   .54,  3.54, -6.46,-12.46, -5.46,
   3.54,  5.54, -8.46, -7.46, -4.46,   .54,  3.54, -9.46,  3.54,  1.54,
   -.46,-12.46, -5.46,-11.46,  4.54,  5.54, -3.46, -9.46,  6.54, -6.46,
 -11.46, -9.46,  9.54,  3.54, -2.46,   .54)

rw <- c(
     4.7859, 6.4711, 5.6384, 4.5445, 5.2512, 5.4320, 6.8870,
     7.3520, 7.5554, 9.4746, 8.8956, 7.8375, 7.1667, 8.0101,
     8.3214, 6.9438, 6.8656, 6.6391, 6.1159, 5.3716, 5.7850,
     5.4549, 5.6714, 5.5356, 4.6291, 5.5859, 4.1332, 5.1065,
     5.1631, 3.7457, 4.3934, 4.1329, 3.5110, 3.6592, 3.4896,
     5.4150, 5.9904, 6.8264, 6.5395, 5.7517, 7.4551, 8.1824,
     7.6346, 6.9486, 6.9439, 7.9076, 8.1756,10.0101,12.5119,
    12.4833,10.8089,12.1630,11.2429,12.4282,11.3600,11.4497,
    10.7938,10.8928, 9.8986,10.4068,10.8555,10.6883,12.9060,
    12.2219,12.2337,13.4497,13.5710,12.6664,11.5258,11.3055,
    10.5535, 9.8174, 9.9523, 9.1888,10.7718,11.3440,10.5945,
     9.0351, 8.0824, 6.9112, 5.7016, 7.4199, 7.0175, 6.8965,
     6.9900, 8.8312, 9.5003, 8.4917, 6.7446, 6.0545, 7.4966,
     7.7099, 6.3240, 5.4536, 4.0820, 5.0447, 3.9409, 3.7853,
     4.2093, 6.1621, 5.8295, 5.7880, 7.9112, 7.2202, 7.4522,
     6.8876, 6.9819, 6.7677, 7.6116, 8.3178, 8.8387, 9.0486,
     8.6076, 7.6684, 8.1584, 7.4693, 5.6946, 4.8324, 6.6357,
     6.3943, 5.7154, 6.5953, 6.1226, 6.0727, 6.9636, 7.6535,
     6.6300, 7.9196, 6.7333, 6.4188, 7.1246, 7.4150, 7.7606,
     7.5092, 8.5348, 8.4524, 7.4116, 6.3278, 7.4545, 7.9191,
     8.9072, 8.9628, 9.1451, 9.7785,10.0496, 9.8638,10.9912,
    11.7322,11.7365,13.5885,13.2428,15.0027,15.2795,14.6536,
    13.8078,13.4062,12.4716,13.3632,13.6777,13.4282,15.5994,
    16.6311,16.2653,16.0370,15.5422,14.7332,14.2666,13.9666,
    15.0999,14.1622,14.4308,12.3841,13.1813,12.8510,12.5237,
    11.9447,11.5721,12.6238,12.2099,13.1526,14.0196,14.3557,
    11.7496,12.4243,12.8397,11.3166,12.2531,14.1968,15.2011,
    14.5817,13.9441,14.6993,15.1034,15.4357,16.3698,16.1636,
    16.5665,15.5200,15.6034,16.5936)

sales <- c(
 154, 96, 73, 49, 36, 59, 95,169,210,278,298,245,
 200,118, 90, 79, 78, 91,167,169,289,347,375,203,
 223,104,107, 85, 75, 99,135,211,335,460,488,326,
 346,261,224,141,148,145,223,272,445,560,612,467,
 518,404,300,210,196,186,247,343,464,680,711,610,
 613,392,273,322,189,257,324,404,677,858,895,664,
 628,308,324,248,272)

star <- c(
 25,28,31,32,33,33,32,31,28,25,22,18,14,10, 7, 4, 2, 0, 0, 0,
  2, 4, 8,11,15,19,23,26,29,32,33,34,33,32,30,27,24,20,17,13,
 10, 7, 5, 3, 3, 3, 4, 5, 7,10,13,16,19,22,24,26,27,28,29,28,
 27,25,24,21,19,17,15,13,12,11,11,10,10,11,12,12,13,14,15,16,
 17,18,19,19,19,19,20,20,20,20,20,20,20,20,21,20,20,20,20,19,
 18,17,16,15,13,12,11,10, 9, 9,10,10,11,12,14,16,19,21,24,25,
 27,28,29,29,28,27,25,23,20,17,14,11, 8, 5, 4, 2, 2, 2, 4, 6,
  9,12,16,19,23,27,30,32,33,34,33,32,30,27,24,20,16,12, 9, 5,
  3, 1, 0, 0, 1, 3, 6, 9,13,17,21,24,27,30,32,33,33,32,31,28,
 25,22,19,15,12, 9, 7, 5, 4, 4, 5, 5, 7, 9,12,14,17,20,22,24,
 25,26,27,27,26,25,24,22,20,18,17,15,14,13,13,12,12,12,13,13,
 13,14,14,15,15,16,17,17,17,17,18,18,19,19,20,20,21,21,22,22,
 22,22,22,21,20,19,17,16,14,12,11, 9, 8, 7, 8, 8, 9,10,12,14,
 17,20,23,25,27,29,30,30,30,29,27,25,22,19,16,12, 9, 6, 4, 2,
  1, 1, 2, 4, 7,10,14,17,21,25,29,31,33,34,34,33,31,29,26,22,
 19,14,11, 7, 4, 2, 1, 0, 1, 2, 5, 7,11,15,19,22,25,28,30,32,
 32,32,31,29,26,23,21,17,14,11, 9, 7, 6, 5, 6, 6, 7, 9,11,13,
 15,18,20,22,23,24,25,25,25,24,24,22,21,19,18,17,16,15,15,14,
 14,14,14,14,14,14,14,14,14,14,15,15,15,15,16,16,17,18,19,20,
 21,22,23,23,24,24,24,23,22,21,19,17,15,13,11, 9, 7, 6, 6, 6,
  7, 8,10,12,15,18,22,24,27,29,31,31,31,31,29,27,24,21,18,14,
 10, 7, 5, 2, 1, 0, 1, 2, 5, 8,12,15,19,23,27,30,32,34,34,34,
 32,30,28,24,20,16,13, 9, 6, 3, 2, 1, 1, 2, 4, 6, 9,13,17,20,
 23,26,28,30,31,31,31,29,27,24,22,19,16,13,11, 9, 8, 7, 7, 7,
  8, 9,11,12,14,16,18,20,21,22,23,23,23,23,23,22,21,20,19,18,
 18,17,17,16,16,16,16,15,15,15,14,14,13,13,13,13,13,13,14,14,
 15,16,18,19,21,22,24,24,25,26,26,25,24,23,21,19,16,14,12, 9,
  7, 5, 5, 4, 5, 6, 8,10,13,16,20,23,26,29,31,32,32,32,31,29,
 26,23,20,16,12, 8, 6, 3, 1, 0, 0, 1, 3, 6,10,13,17,21,25,28,
 31,33,34,34,33,31,29,26,22,18,15,11, 8, 5, 3, 2, 2, 2, 4, 5)

waves <- c(
   -.103, .183, .241, .281, .298, .065, .102, .143,-.029,-.083,
   -.315,-.329,-.258,-.255, .001, .204, .318, .595, .585, .524,
    .315, .059, .005,-.191,-.545,-.693,-.774,-.760,-.407,-.285,
   -.002, .251, .318, .534, .460, .247, .190,-.053,-.012,-.053,
   -.241,-.191,-.248,-.295,-.171,-.154, .018, .217, .183, .298,
    .227, .112,-.096,-.359,-.447,-.481,-.578,-.417,-.339,-.225,
    .197, .311, .423, .402, .143, .038,-.002,-.117,-.083,-.161,
   -.268,-.127,-.123, .180, .375, .308, .399, .173,-.093,-.154,
   -.393,-.305,-.208,-.046, .406, .621, .389, .220,-.204,-.299,
   -.356,-.390,-.315,-.326,-.154,-.120,-.174, .089, .278, .328,
    .440, .352, .096, .059,-.218,-.299,-.346,-.410,-.285,-.275,
   -.194, .150, .268, .402, .547, .416, .345, .042,-.069,-.093,
   -.322,-.373,-.329,-.258,-.069, .055, .204, .251, .109, .089,
    .166, .045, .072,-.063,-.245,-.191,-.157,-.009, .069, .005,
    .069, .005,-.103, .015,-.100,-.019, .140, .001, .008,-.100,
   -.262,-.302,-.403,-.272,-.059,-.049, .072,-.016,-.171,-.157,
   -.258,-.285,-.133, .008, .055, .230, .234, .416, .487, .402,
    .305, .022,-.150,-.315,-.481,-.504,-.760,-.730,-.521,-.386,
   -.019, .160, .214, .426, .571, .726, .928, .830, .601, .197,
   -.191,-.164,-.494,-.706,-.642,-.471,-.181, .008, .177, .507,
    .591, .618, .460, .284, .271, .112,-.120,-.177,-.187,-.326,
   -.225,-.265,-.140,-.022,-.127, .112, .251, .301, .409, .237,
    .109, .005,-.167,-.043,-.073,-.133, .048,-.063, .011, .160,
    .136, .237, .268, .112, .055,-.106,-.181,-.225,-.356,-.161,
   -.063,-.076, .092, .109, .123, .288, .227, .349, .500, .180,
   -.039,-.275,-.450,-.467,-.619,-.457,-.292,-.137, .129, .399,
    .426, .493, .268, .281, .160,-.127,-.204,-.433,-.605,-.488,
   -.430,-.059, .244, .551, .884, .638, .146,-.029,-.332,-.352,
   -.346,-.572,-.339,-.245,-.332,-.238,-.059, .217, .244, .086,
   -.022,-.154,-.309,-.322,-.444,-.275,-.198,-.241,-.016, .079,
    .099, .284, .065, .055, .062,-.090,-.026,-.049,-.127, .109,
    .001,-.009,-.056,-.127, .008, .059, .075, .187, .106, .005)

wolfer <- c(
 80.9,  83.4,  47.7,  47.8,  30.7,  12.2,   9.6,  10.2,  32.4,  47.6,
 54.0,  62.9,  85.9,  61.2,  45.1,  36.4,  20.9,  11.4,  37.8,  69.8,
106.1, 100.8,  81.6,  66.5,  34.8,  30.6,   7.0,  19.8,  92.5, 154.4,
125.9,  84.8,  68.1,  38.5,  22.8,  10.2,  24.1,  82.9, 132.0, 130.9,
118.1,  89.9,  66.6,  60.0,  46.9,  41.0,  21.3,  16.0,   6.4,   4.1,
  6.8,  14.5,  34.0,  45.0,  43.1,  47.5,  42.2,  28.1,  10.1,   8.1,
  2.5,   0.0,   1.4,   5.0,  12.2,  13.9,  35.4,  45.8,  41.1,  30.4,
 23.9,  15.7,   6.6,   4.0,   1.8,   8.5,  16.6,  36.3,  49.7,  62.5,
 67.0,  71.0,  47.8,  27.5,   8.5,  13.2,  56.9, 121.5, 138.3, 103.2,
 85.8,  63.2,  36.8,  24.2,  10.7,  15.0,  40.1,  61.5,  98.5, 124.3,
 95.9,  66.5,  64.5,  54.2,  39.0,  20.6,   6.7,   4.3,  22.8,  54.8,
 93.8,  95.7,  77.2,  59.1,  44.0,  47.0,  30.5,  16.3,   7.3,  37.3,
 73.9, 139.1, 111.2, 101.7,  66.3,  44.7,  17.1,  11.3,  12.3,   3.4,
  6.0,  32.3,  54.3,  59.7,  63.7,  63.5,  52.2,  25.4,  13.1,   6.8,
  6.3,   7.1,  35.6,  73.0,  84.9,  78.0,  64.0,  41.8,  26.2,  26.7,
 12.1,   9.5,   2.7,   5.0,  24.4,  42.0,  63.5,  53.8,  62.0,  48.5,
 43.9,  18.6,   5.7,   3.6,   1.4,   9.6,  47.4,  57.1, 103.9,  80.6,
 63.6,  37.6,  26.1,  14.2,   5.8,  16.7,  44.3,  63.9,  69.0,  77.8,
 65.0,  35.7,  21.2,  11.1,   5.7,   8.7,  36.1,  80.4, 114.4, 109.6,
 88.8,  67.8,  47.5,  30.6,  16.3,   9.6,  33.2,  92.6, 151.6, 136.3,
134.7,  83.9,  69.4,  31.5,  13.9,   4.4,  38.0, 141.7, 190.2, 184.8,
159.0, 112.3,  53.9,  37.5,  27.9)

  return(list(air=air,artif=artif,bev=bev,buffsnow=buffsnow,
              calfem=calfem,cos4=cos4,cradfq=cradfq,
              cradfqsun=cradfqsun,eriel=eriel,gasfurnin=gasfurnin,
              gasfurnout=gasfurnout,
              lh1=lh1,lh2=lh2,lynx=lynx,mlco2=mlco2,normwn=normwn,
              nyct=nyct,nycb=nycb,raineast=raineast,rw=rw,
              sales=sales,star=star,waves=waves,wolfer=wolfer))
}
#
#
#
windowf <- function(rho, R0, Q, ioptw, M, n, alpha=.05)
#---------------------------------------------------------------------
#
#   Function to find  nonparametric spectral estimator 
#   at Q frequencies in [0,1]
#   using window number ioptw, scale parameter M. 
#   The sample size n is also input. 
#
#   The spectral estimator and the constant c used in confidence
#   intervals are returned in a list.
#
#---------------------------------------------------------------------
{
  
  if(M >= n || Q < n || ioptw < 1 || ioptw > 8 
     || (ioptw <= 5 && length(rho) < M)
     || (ioptw > 5 && length(rho) < n - 1)
     || alpha <= 0. || alpha >= 1.)
    stop("Illegal Input to window()")
  
  z          <- rep(0, Q)
  if(ioptw <= 5) z[1:(M+1)] <- c(R0, R0 * rho[1:M])
  if(ioptw > 5)  z[1:n] <- c(R0, R0 * rho[1:(n-1)])
  
  z1 <- rep(0,Q)
  
  if(ioptw == 1) z1[1:(M+1)] <- rep(1, M+1)
  
  else if(ioptw == 2)  z1[1:(M+1)] <- 1 - (c(0:M) / M)
  
  else if(ioptw == 3)  z1[1:(M+1)] <- .54+.46*cos(pi*c(0:M)/M)
  
  else if(ioptw == 4) {
    u  <- c(1:M)/M
    u1 <- u[u<=.5]
    u2 <- u[u>.5]
    z1[1:(M+1)]  <- c(1,1-6*u1^2+6*u1^3,2*(1-u2)^3)
  }
  
  else if(ioptw == 5) {
    u   <- c(0:M) / M
    piu <- pi * u
    z1[1:(M+1)]  <- (1 - u) * cos(piu) + sin(piu) / pi
  }
  
  else if(ioptw == 6) {
    piu <- pi * (c(1:(n-1)) / M)
    z1[1:n] <- c(1,sin(piu)/piu)
  }
  
  else if(ioptw == 7) {
    piu <- pi * (c(1:(n-1)) / M)
    z1[1:n] <- c(1,3*((sin(piu)/piu) - cos(piu)) / (piu*piu))
  }
  
  else if(ioptw == 8) {
    u <- c(0:(n-1))/M
    z1[1:n] <- 1 /(1+u^4)
  }
  
  z <- c(z1[1] * z[1], 2 * z1[2:Q] * z[2:Q])
  
  z <- Re(fft(z))
  z <- z[1:((Q / 2) + 1)]
  
  ilam2 <- c(2,2./3.,.795,.539,.586,1,1.2,1.66)
  fac <- qnorm(1-(alpha/2)) * sqrt(M*ilam2[ioptw]/n)
  
  return(list(f=z,c=fac))
  
}

windsp3 <- function(x,ioptw=4,Q=length(x),M1=10,M2=20,M3=40,   #
                    main="Log Standardized Window Spectral Estimates")#
  #--------------------------------------------------------
#
#  windsp3: plot log standardized window spectral estimator
#           of x using window specified by ioptw and three
#           truncation points M1 <= M2 <= M3 at the fequencies
#           (j-1)/Q, j=1,...,[Q/2]+1.
#
#---------------------------------------------------------
{
  
  if(ioptw < 1 || ioptw >5) {
    print("ioptw must be between 1 and 5")
    return()}
  
  
  wnds <- c("Truncated periodogram","Bartlett","Tukey","Parzen",
            "Bohman")
  
  n <- length(x)
  
  rho=corr(x,M3)
  
  f1 <- windowf(rho$corr,rho$var,Q,ioptw,M1,n)$f
  f2 <- windowf(rho$corr,rho$var,Q,ioptw,M2,n)$f
  f3 <- windowf(rho$corr,rho$var,Q,ioptw,M3,n)$f
  
  plotsp(f3,Q,rho$var,main=main)
  lines(freqs(Q),log(f2/rho$var),lty=2)
  lines(freqs(Q),log(f1/rho$var),lty=3)
  
#  text(0,5.5,paste(wnds[ioptw]," window"),adj=0)
#  text(0,5,paste("M1 = ",M1," M2 = ",M2," M3 = ",M3),adj=0)
  legend("topleft",lty = c(1,2,3), bty="n",ncol=1,cex=0.85,title=paste(wnds[ioptw],"window"),
         legend = c(as.expression(bquote(M[1]==.(M1))),
                    as.expression(bquote(M[2]==.(M2))),
                    as.expression(bquote(M[3]==.(M3)))))
#  
  return(invisible())
}