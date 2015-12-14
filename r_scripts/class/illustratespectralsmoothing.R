#
#  AR spectrum: 
#
par(mfrow = c(2, 1))
alpha <- c(0.1, 0, 0, -0.8)
n     <- 256
Q     <- 256
x     <- ardt(alpha = alpha, n = n)$x
div   <- arcorr(alpha = alpha)
f     <- arsp(alpha = alpha, Q = Q)
fhat  <- perdgm(x)
plotsp(f = f, n = n, div = div$var, main = "AR Spectrum")
lines(freqs(Q), log(stdf(f = fhat, fac = div$var, a = exp(-6), b = exp(6))))
windsp3(x = x)
lines(freqs(Q), log(stdf(f = fhat, fac = div$var, a = exp(-6), b = exp(6))))
#
#  MA spectrum:
#
par(mfrow = c(2, 1))
beta <- c(0.1, 0, 0, -0.8)
n    <- 256
Q    <- 256
x    <- madt(beta = beta, n = n)
div  <- macorr(beta = beta)
f    <- masp(beta = beta, Q = Q)
fhat <- perdgm(x)
plotsp(f = f, n = n, div = div$var, main = "MA Spectrum")
lines(freqs(Q), log(stdf(f = fhat, fac = div$var, a = exp(-6), b = exp(6))))
windsp3(x = x)
lines(freqs(Q), log(stdf(f = fhat, fac = div$var, a = exp(-6), b = exp(6))))
#
#  ARMA spectrum:
#
par(mfrow = c(2, 1))
alpha <- c(0, 0, -0.9)
beta  <- c(0, 0, 0, 0.8)
n     <- 256
Q     <- 256
x     <- armadt(alpha = alpha, beta = beta, n = n)$x
div   <- armacorr(alpha = alpha, beta = beta)
f     <- armasp(alpha = alpha, beta = beta, Q = Q)
fhat  <- perdgm(x)
plotsp(f = f, n = n, div = div$var, main = "ARMA Spectrum")
lines(freqs(n), log(stdf(f = fhat, fac = div$var, a = exp(-6), b = exp(6))))
windsp3(x = x)
lines(freqs(n), log(stdf(f = fhat, fac = div$var, a = exp(-6), b = exp(6))))