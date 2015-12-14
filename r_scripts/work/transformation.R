# install.packages("zoo")
library(zoo)
source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')

# i'th degree polynomial
y = 2 + 5*t+0.5*t^4
plotl(diff(diff(diff(diff(y, 1), 1))))


# Transformation of air TS
plotl <- function(x, ylab="") plot(x, type="l", ylab=ylab)
air = datasets()$air
t = 1:length(air)
par(mfrow=c(2,3), mar = c(2,4,2,1) + 0.1)
plotl(air, "air")
plotl(log(air), "log air") # to combat heteroscedasticity
plotl(diff(air,1), "diff 1 air") # hetero
plotl(diff(log(air), 1), "diff 1 of log air")
plotl(lm(log(air)~t)$residuals, "residuals of lm(log(air)~t)") # remove linear trend
plotl(diff(diff(log(air), 1), 12), "Monthly diff of diff 1 of log air") # remove seasonality


# Seasonal Variability
plotl(submns(air, d = 12)$x)
plotl(divsds(submns(air, d = 12)$x, d = 12)$x)

stdts <- function(data, d){
    xbar = vector(); ss = vector();
    for (start in 1:d){
        ints = seq(start, length(data), by=d)
        x = data[ints] 
        m = mean(x); xbar = c(xbar, m);
        s = sd(x); ss = c(ss, s);
        data[ints] = (x - m) / s
    }   
    return (list(x=data, xbar=xbar, sds=ss))
}
# Test
plot(stdts(air, 12)$x, type="l")
plot(stdts(air, 12)$sds, type="l")
plot(stdts(air, 12)$xbar, type="l")

# Same
plot(submns(air, 12)$xbar, type="l")
plot(divsds(air, 12)$sds, type="l")

# Season Adjusted, trend removed
plot(diff(stdts(air, 12)$x, 1), type="l")
wntest(diff(stdts(air, 12)$x, 1)) # Still something going on


# Moving Average Smoother
mav <- function(x,n=5){na.omit(filter(x,rep(1/n,n), sides=2))}
gkg = read.table('gkg_10032.csv')$V1
par(mfrow=c(3,1), mar = c(2,4,2,1) + 0.1)
plotl(gkg)
plotl(mav(gkg, n = 5))
plotl(rollmean(gkg, 5))

### Inverse Differenciating Forecast ###
plotl(extend(air, h=24, d1=1, d2=12))
points(air)

### Simple Exponential Smoothing and Forecasting ###
expsm(rw)$alpha