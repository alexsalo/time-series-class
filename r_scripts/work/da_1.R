source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/sawtooth.R')
source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')

######################################################
############# Data Analysis: Assignment 1#############
######################################################

# Read Data
ts1 = read.table('series1.dat')$V1
ts2 = read.table('series2.dat')$V1
ts3 = read.table('series3.dat')$V1
ts4 = read.table('series4.dat')$V1
ts5 = read.table('series5.dat')$V1
ts6 = read.table('series6.dat')$V1
ts7 = read.table('series7.dat')$V1
gkg = read.table('gkg_10032.csv')$V1

# TS 1
descplot(ts1, 40)
wntest(ts1, 40)

# TS 2
descplot(ts2, 40)
wntest(ts2, 40)

# TS 3
ts3 = ts3 - ts1[1:240]
descplot(ts3, 40)
plotsp(perdgm(ts3), length(ts3), sd(perdgm(ts3)))
ts3.freqs = round(Mod(fft(ts3))[1:(length(ts3)/2+1)] * 2 / length(ts3), 6)
ts3.freqs[ts3.freqs < 1] = 0

ts3.freqs[ts3.freqs > 0] # amps
length(ts3) / (-1 + which(ts3.freqs > 0)) #freqs

par(mfrow=c(2,1),mar = c(4,5,1.5,1) + 0.1)
t = 1:240
plot(1.2*cos(t*2*pi/80) + 3.8*cos(t*2*pi/12) + rnorm(240), type="l")
plot(ts3, type="l")

# TS 4
ts4 = ts4[1:240] - ts1[1:240]
descplot(ts4, 40)
ts4.freqs = round(Mod(fft(ts4))[1:(length(ts4)/2+1)] * 2 / length(ts4), 6)
ts3.freqs[ts4.freqs < 1] = 0
ts4.freqs[ts4.freqs > 0] # amps
length(ts4) / (-1 + which(ts4.freqs > 0)) #freqs

# TS5
descplot(ts5, 40)
wntest(ts5)
AR <- function(n=250, coef=c(0, 0, 0)){
	x = c(1,1,1)
	noise = rnorm(n)
	for (i in 4:n){
		x = c(x, coef[1]*x[i-1] + coef[2]*x[i-2] + coef[3]*x[i-3]+ noise[i-3])
	}
	return(x[3:(n+3])
}	

ar5 = AR(n=250, coef=c(0, 0.73, 0))
descplot(ar5, 40)

descplot(arima.sim(model=list(ar=c(0.15, -0.73)), n=250), 40)

# TS6
descplot(ts6, 40)
MA <- function(n=250, coef=c(0, 0, 0)){
	noise = rnorm(n + 3)
	x = coef[1] * (noise[4:(n+3)] - noise[3:(n+2)]) + 
		coef[2] * (noise[4:(n+3)] - noise[2:(n+1)]) + 
		coef[3] * (noise[4:(n+3)] - noise[1:n])
	return(x)
}
ma6 = MA(n=250, coef=c(0,1.73, 0))
descplot(ma6, 40)

descplot(arima.sim(model=list(ma=c(0.15, -0.73)), n=250), 40)

# TS7
descplot(ts7, 40)
descplot(arima.sim(model=list(ar=c(-0.3, 0.7, 0.2), ma=c(0, 0.2)), n=250), 40)

ARMA <- function(n=250, arc=c(0, 0, 0), mac=c(0, 0, 0)){
	x = c(1,1,1)
	noise = rnorm(n + 3)
	for (i in 4:(n+3)){
		x = c(x, arc[1]*x[i-1] + arc[2]*x[i-2] + arc[3]*x[i-3]+ noise[i])
	}
	noise = mac[1] * (noise[4:(n+3)] - noise[3:(n+2)]) + 
	mac[2] * (noise[4:(n+3)] - noise[2:(n+1)]) + 
	mac[3] * (noise[4:(n+3)] - noise[1:n])
	
	return(x[4:(n+3)] + noise)
}
arma7 = ARMA(arc=c(-0.3, 0.7, 0.15), mac=c(0, 0.6, 0))
descplot(arma7, 40)

# GKG
descplot(gkg, 40)
plotl(rollmean(gkg, 5))
descplot(rollmean(gkg, 5), 40)
wntest(rollmean(gkg, 5), 40)

gkg = rollmean(gkg, 5)
gkg.freqs = round(Mod(fft(gkg))[1:(length(gkg)/2+1)] * 2 / length(gkg), 6)
gkg.freqs[gkg.freqs < 0.3] = 0
gkg.freqs[gkg.freqs > 0] # amps
length(gkg) / (-1 + which(gkg.freqs > 0)) #freqs