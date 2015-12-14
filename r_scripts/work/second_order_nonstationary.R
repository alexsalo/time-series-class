source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')

x = rw(200)
t = 1:200
par(mfrow=c(3,1))
plotl(x)
plotl(lm(x ~ t)$residuals)
plotl(diff(x,1)) # undo random walk, leave noise only

# Non Second Order Stationary
N = 10
data = vector()
for (i in 1:N){ # Creating the Ensemble
	data = rbind(data, rw(200))
}	
max = max(data)
min = min(data)
R = 0.3; G = 0.4; B = 0.9;
color = rgb(R, G, B, 0.6)
plot(data[1,1:200], ylim=c(min, max), type="l", col=color)
for (i in 2:N){
	greens = seq(0, 0.5, length.out=N-1)
	color = rgb(R, G+greens[i], B, 0.4)
	lines(data[i,1:200], type="l", col=color)
}	