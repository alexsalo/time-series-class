install.packages("gridExtra")
install.packages("smoothmest")
install.packages("ggplot2")
library(ggplot2)
library(gridExtra)
library(smoothmest)
source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')

##########################################
############### HOMEWORK 1 ###############
##########################################

# Plot 2x2 TS from DF
plot_wn <- function(df){
	ggplot(data=df, aes(x=time, y=value)) + 
		geom_line(lwd = 0.5) + labs(x = "Time", y = "Value")
}

gen_df <- function(dist_fun, plotname){
	plot.new()
	ns = c(50, 100, 250, 500)
	df = data.frame(time = seq(ns[4]), value = wn(123, ns[4], dist_fun)) #dist_fun(ns[4]))	
	do.call(grid.arrange, c(lapply(ns, function(x) plot_wn(df[1:x,])), main=plotname))
}

gen_df(1, "White Noise") # C1.1
gen_df(2, "Uniform") # C1.2
gen_df(7, "Log-Normal") # C1.7
gen_df(8, "Double Exponential") # C1.8

### C1.9
N = 501
thetas = seq(0, 2*pi, length.out = N)

plot_thetas <- function(data, ylabel){
	df = data.frame(time = thetas, values = data)
	ggplot(data=df, aes(x=time, y=values)) + geom_line() +
	labs(x = "t", y = ylabel) + ylim(-2, 2) +
	scale_x_continuous(breaks=c(0,pi,2*pi), labels=c("0", pic, pic2))
}

grid.arrange(plot_thetas(cos(thetas), "cos(x)"),
			 plot_thetas(cos(2 * thetas), "cos(2x)"),
			 plot_thetas(cos(thetas / 2), "cos(x/2)"),
			 plot_thetas(2 * cos(thetas), "2*cos(x)"),
			 plot_thetas(2 * cos(2 * thetas), "2*cos(2x)"),
			 plot_thetas(2 * cos(thetas / 2), "2*cos(x/2)"),
 			 ncol=3)		
# Multiply x by const -> increase (dec) frequency; cos * const -> inc(dec) amplitude.

### C1.10
t = seq(0, 100, length.out=1000)
pic = expression(pi); pic2 = expression(2*pi);
df = data.frame(time = t, sine = sin(t*pi/10), cosine = cos(t*pi/10), 
	sincos = sin(t*pi/10) + cos(t*pi/10))
ggplot(data=df) + geom_line(aes(x=time, y=sine, colour="sin(x)")) + 
				  geom_line(aes(x=time, y=cosine, colour="cos(x)")) + 
				  geom_line(aes(x=time, y=sincos, colour="sin(x) + cos(x)")) #+
		#scale_x_continuous(breaks=c(0,pi,2*pi), labels=c("0", pic, pic2))

### T1.4(b)
hist(rnorm(100)/rnorm(100))
plot(atan(rnorm(100)/rnorm(100))) 	

##########################################
############### HOMEWORK 2 ###############
##########################################

# TimesLab.R
par(mfrow=c(3,2))
plot(wn(100), type = "l", main = "White Noise")
plot(datasets()$air, type = "l", main = "Air Trafic Data")

plotcorr(corr(wn(100), 30)$corr, "Correlogramm of White Noise")
plotcorr(corr(datasets()$air, 30)$corr, "Correlogramm of Air Trafic Data")

plotcorr(parcorr(wn(100), 30), "Partial ACF of White Noise")
lines(srvars(wn(100), 30), type="p")
plotcorr(parcorr(datasets()$air, 30), "Partial ACF of Air Trafic Data")
lines(srvars(datasets()$air, 30), type="p")

# R built-in functions
par(mfrow=c(2,2))
acf(wn(100), 30, "correlation")
acf(datasets()$air, 30, "correlation")
acf(wn(100), 30, "partial")
acf(datasets()$air, 30, "partial")

# Descriptive
descplot(datasets()$air,30)

# Periodogram Sandbox
z = round(Mod(fft(datasets()$air)[1:101]),2)
plot(z, type="line")

# Alex Salo ACF Plots
acfplots <- function(time, data, numlags, name="Correlogramm"){
	par(mfrow=c(3,1),mar = c(2,5,1.5,1) + 0.1)
	plot(time, data, type = "l", main = name)
	plotcorr(corr(data, numlags)$corr, paste("Correlogramm of", name, sep=" "))
	plotcorr(parcorr(data, numlags), paste("Partial ACF of", name, sep=" "))
}

acfplots <- function(time, data, numlags, name="Correlogramm"){
	par(mfrow=c(2,1),mar = c(2,5,1.5,1) + 0.1)
	plot(time, data, type = "l", main = name)
	plot(corr(x, 40)$corr, type="l", main  = paste("Correlogramm of", name, sep=" "))
}

### C2.3
random_walk <- function(n=200){
	a = wn(n,"unif")
	b = rep(-1, n)
	b[a > 0.5] = 1
	return(cumsum(b))	
}
rw = random_walk(200)
acfplots(1:200, rw, 40, "Random Walk") 
# Since it's long term memory, we expect the auto-correlation to be high at the small lags
# and evade with the increase of the lag size.

### C2.6
stooth <- function(ncycles = 5, amp = 10, plot.it = as.logical(T))
{
	n <- 4*amp*ncycles
	t <- 1:n
	x <- rep(0, n)
	a2 <- 2*amp
	
	x[1:amp] <- t[1:amp]
	x[(n-amp+1):n] <- -amp + t[1:amp]

	for(cycle in 1:(2*ncycles-1)) {
		cc <- 2*cycle - 1
		dd <- 2*cycle + 1
		x[(cc*amp + 1):(dd*amp)] <- ((-1)^(cycle+1))*amp + ((-1)^cycle)*t[1:a2]
	}
	
	if(plot.it)
		plot(x, type = "l", xlab = "t", ylab = expression(x[t]), main = "Sawtooth Function")
	
	return(list(t = t, x = x))
}

sawtooth <- function(ncycles = 5, amp = 10){
	precision = 2 * amp
	t = seq(0, 2 * amp, length.out = precision)
	cycle = c(t[1:length(t)-1], rev(t)[1:length(t)-1]) - amp
	cycle = c(cycle[(precision / 2): length(cycle)], cycle[0: (precision / 2)])
	cycles = rep(cycle[2:length(cycle)], ncycles)	
	return(cycles)
} 
st = stooth(12, 2)$x
acfplots(1:length(st), st, 40, 'Saw Tooth')
descplot(st, 40)

which(round(perdgm(st), 3) > 0)
plotsp(perdgm(st), length(st), sd(perdgm(st)))


### C2.7
t = seq(0, 100)
x = sin(t * pi / 2.5)
acfplots(t, x, 40, "Sine Wave")

plot(corr(x, 40)$corr, type="l")
# So correlogram of a Sine is also a Sine. R(0) = 1. 

### C2.8
plot_superimposed_sines <- function(n, period1, period2){
	t = seq(0, n)
	x1 = sin(t * pi / (period1 / 2))
	x2 = sin(t * pi / (period2 / 2))
	x3 = x1 + x2
	plot(t, x3, type="l", col="red", lwd="2")
	lines(t,x1,type="l", col="green")
	lines(t,x2,type="l", col="blue")
	legend(0, -1.2, c("sin(p1)", "sin(p2)", "sin(p1) + sin(p2)"), lty=c(1,1,1),
		lwd=c(2.5,2.5,2.5),col=c("green","blue","red"))
	
	dot(x1, x2)
}
plot_superimposed_sines(100, 10, 20)
plot_superimposed_sines(100, 10, 10)

### C2.12
primes = c(19991, 19993, 199999)
num_sample = c(100, 1000, 1500, 2000, 5000)
timecov <- function(n = 500, M = 40){
	x <- rnorm(n)
	t.corr <- system.time(corr(x, M))[3]
	t.corr1 <- system.time(corr1(x, M))[3]
	return(c(unname(t.corr), round(unname(t.corr1), 2)))
}

for (prime in primes){
	for (M in num_sample){
		time = timecov(primes, M)
		print(paste("Prime: ", prime, " M: ", M, ", 2FFTs: ", time[1], ", Conv: ", time[2], sep=""))
	}
}

system.time(corr(19999999, 5000))

### C2.14
plot_three_sine <- function(n=200, noise_variance=0,
		amplitudes=c(10, 3, 1), frequencies=c(100, 10, 4)){
	t = seq(1, n)	
	result = vector()
	for (i in 1:3){
		result = rbind(result, amplitudes[i] * cos(2*pi*t / frequencies[i]))
	}
	sum = colSums(result)
	if (noise_variance != 0)
		sum = sum + rnorm(n, 0, noise_variance)
	plot(t, sum, type="l")	
	return (sum)
}

plot_std_periodogram <- function(ts){
	N = length(ts)
	variance = sum((ts - mean(ts))^2) / N
	zk_sq = Mod(fft(ts)^2)[1:(N/2+1)]
	Ck_sq = zk_sq / N
	# Standartize and truncate
	Ck_sq = Ck_sq/variance
	Ck_sq[Ck_sq < exp(-6) ] = exp(-6)
	Ck_sq[Ck_sq > exp(6) ] = exp(6)
	k = 1:(N/2+1)
	omega_k = (k - 1) / N
	plot(omega_k, log(Ck_sq), type="l", ylim = c(-6, 6), 
		xlab = expression(omega), ylab = expression(hat(f)*(omega)),
		main = "Natural Log of Standardized Spectrum")
	return(list(frequencies=round(Mod(fft(ts))[1:(N/2+1)] * 2 / N, 6),
 				periodogram=log(Ck_sq)))
}

C2.14 <- function(noise_variance=0,
		amplitudes=c(10, 3, 1), frequencies=c(100, 10, 4)){	
	par(mfrow=c(2,1),mar = c(4,5,1.5,1) + 0.1)
	sine_sum = plot_three_sine(noise_variance=noise_variance,
		amplitudes=amplitudes, frequencies=frequencies)
	periodogram = plot_std_periodogram(sine_sum)
	return(list(sine_sum=sine_sum, frequencies=periodogram$frequencies, 
			    periodogram=round(periodogram$periodogram, 2)))
}
c214 = C2.14(0, c(10, 3, 1))

# reconstruct sine waves
N = length(c214$sine_sum)
amps = c214$frequencies[c214$frequencies > 0]
freqs = N / (-1 + which(c214$frequencies > 0))
plot_three_sine(amplitudes=amps, frequencies=freqs)

### C2.17
variance = c(0.01, 0.25, 4)
C2.14(noise_variance=variance[1])
C2.14(noise_variance=variance[2])
C2.14(noise_variance=variance[3])


# GKG 10032
gkg = read.table('gkg_10032.csv')$V1
psych::describe(gkg)
acfplots(seq(0, length(gkg) - 1), gkg, 40, "GKG")

### T2.1
t = seq(1,50)
x = 5 + 1*t
res = corr(x, 50)
steps = res$corr[2:length(res$corr)] - res$corr[1:(-1+length(res$corr))]
plot(steps, type="l", xlab="lag", ylab="Step size = corr_(v+1) - corr_v", 
	main="Steps between correlogram of a line")
acfplots(t,x,20,"l")
plot(corr2(x[1:70], x[31:100], 30)$cross.corr[31:60])

x = 3*(-t^4 + 2 * 50 * t^3 + 2500 * t^2)/(50^3 - 50)
plot(x, main="-t^4 + 2 * n * t^3 + n^2 * t^2")