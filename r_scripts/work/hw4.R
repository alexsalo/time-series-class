source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')
install.packages('knitr', dependencies = TRUE)
library(knitr)

### C4.2
gen_wn_corr <- function(n, nreps=20, m=25, sig2=1)
{
	data = vector()
	for (i in 1:nreps){ # Creating the Ensemble
		data = rbind(data, corr(rnorm(n,0,sig2), m)$corr)
	}	
	max = max(data); max = 0.4
	min = min(data); min = -0.4
	R = 0.3; G = 0.4; B = 0.9;
	color = rgb(R, G, B, 0.6)
	plot(data[1,1:m], ylim=c(min, max), type="l", xlab = "v", col=color, 
		ylab = expression(hat(rho[v])), main = "Correlogram")
	for (i in 2:nreps){
		greens = seq(0, 0.5, length.out=nreps)
		color = rgb(R, G + greens[i], B, 0.6)
		lines(data[i,1:m], type="l", col=color)
	}
	abline(h=0.0, col = 'red', lwd=1.5)
	means = colMeans(data)
	lines(means, type="l", lwd=2)
	return(means)
}

cors.mean = gen_wn_corr(75, nreps=20, m=25, sig2=1)
mean(cors.mean)

par(mfrow=c(3,1),mar = c(2,2,2,2))
for (n in c(30, 250, 2500)) {
	cors.mean = gen_wn_corr(n, nreps=20, m=25, sig2=1)
}


### C4.3
gen_wn_perdgm <- function(n, nreps=20, sig2=1)
{
	lognorm = function(x){log(stdf(x, sig2, exp(-6), exp(6)))}
	data = vector()
	for (i in 1:nreps){ # Creating the Ensemble
		data = rbind(data, perdgm(rnorm(n,0,sig2)))
	}	
	R = 0.3; G = 0.4; B = 0.9;
	color = rgb(R, G, B, 0.6)
	m = (n/2)+1
	plot(freqs(n), lognorm(data[1,1:m]), type="l", col=color,
		xlim = c(0, 0.5), ylim = c(-6, 6), main = "Periodogram",
		xlab = expression(omega), ylab = expression(hat(f)*(omega)))
	for (i in 2:nreps){
		greens = seq(0, 0.5, length.out=nreps)
		color = rgb(R, G + greens[i], B, 0.6)
		lines(freqs(n), lognorm(data[i,1:m]), type="l", col=color)
	}
	abline(h=lognorm(sig2), col = 'red', lwd=1.5)
	means = colMeans(data)
	lines(freqs(n), lognorm(means), type="l", lwd=2)
	return(means)
}
perd.mean = gen_wn_perdgm(n=75, nreps=20, sig2=1)
mean(perd.mean)

par(mfrow=c(3,1),mar = c(2,2,2,2))
for (n in c(30, 75, 500)) {
	perd.mean = gen_wn_perdgm(n, nreps=20, sig2=5)
}


### C4.4

x = wn(n=200);									 descplot(x, m=20)

y = filt(beta=c(0.9),			  beta0=1, x=x); descplot(y, m=20)
y = filt(beta=c(-0.9),			  beta0=1, x=x); descplot(y, m=20)
y = filt(beta=c(0, 0, 0, 0.8),	  beta0=1, x=x); descplot(y, m=20)
y = filt(beta=c(-0.7, -0.1, 0.6), beta0=1, x=x); descplot(y, m=20)


### C4.5
linearTrendDiff=function(a,b,sig2,n){
	lt = b*seq(1,n) + a + rnorm(n,0,sig2)
	d = diff(lt, 1)
	#
	lt.corr=corr(lt,40)$corr; d.corr=corr(d,40)$corr
	lt.var=corr(lt,40)$var;   d.var=corr(d,40)$var
	lt.pgdm=perdgm(lt);		  d.pgdm=perdgm(d)	
	#
	par(mfrow=c(3,2),mar = c(2,2,2,2))
	plot(lt, type="l", main="Linear Trend Time Plot")
	plot(d, type="l", main="First Diff Time Plot")
	abline(h=b, lwd=2, lty=2)
	#
	lim <- qnorm((1+(1-0.05)^(1/40))/2)/sqrt(n)
	plot(lt.corr, ylim = c(-1,1), type="l", main="Linear Trend Correlogram")
	abline(h=lim, lty=2); abline(h=-lim, lty=2)
	plot(d.corr, ylim = c(-1,1), type="l", main="First Diff Correlogram")
	abline(h=lim, lty=2); abline(h=-lim, lty=2)
	#	
	plotsp(lt.pgdm, n, lt.var)
	plotsp(d.pgdm, n-1, d.var) #n-1 for difference
}
linearTrendDiff(1,2,0.25,200)
