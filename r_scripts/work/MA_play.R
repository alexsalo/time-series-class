source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')
options(scipen=999) # to disable scientific notation


plotcorr <- function(z){
  plot(z$corr, ylim = c(-1,1), type = "l",
       xlab = "v", ylab = expression(hat(rho[v])), 
       main = "Correlogram")
}

plotpart <- function(z){
  plot(z$theta, ylim = c(-1, 1), type = "l",
       main = "Partial Correlogram",
       xlab = "v", ylab = expression(hat(theta[v])))
}

trueAR <- function(beta){
  # 1. true parameters
  acor = macorr(beta, 20)
  pcor = mapart(beta, 20)
  sp = masp(beta, 256)
  
  # 2. plot them along with one of the realizations
  par(mfrow=c(2,2))
  tsplot(madt(beta, 200))
  plotcorr(acor)
  plotpart(pcor)
  plotsp(sp, 256, acor$var)
}

trueAR(c(0,0,0,-0.85))
trueAR(c(0.9, 0.25))

# try cumulative SP
sp
plotsp(sp, 256, acor$var)
cumsum(sp)
plot(cumsum(masp(c(0,0,0,-0.85), 256)), type='l')

