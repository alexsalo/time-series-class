library(forecast)

goog = read.csv('C:\\Dropbox\\Baylor\\TimeSeries\\R\\ts_hw\\goog.csv')
names(goog)[1] = "Date"
goog = goog[rev(rownames(goog)),]
plot(goog$Close, type="l")
goog$index = 1:length(goog$Close)

descplot(goog$Close, 100)
par(mfrow=c(1,1))

fit = lm('Close ~ index', goog)
plot(fit)

plot(goog$Close, type="l")
lines(fit$fitted.values, type="l")

goog.defit = goog$Close - unname(fit$fitted.values)
plot(goog.defit, type="l")

goog.detrend = diff(goog.defit, 1)
plot(goog.detrend, type="l")

auto.arima(x=goog.detrend, d=0, D=0, stationary=F, seasonal=F, ic="aicc")
auto.arima(x=goog.detrend, d=0, D=0, stationary=F, seasonal=F, ic="bic")

goog.fit = arima(x=goog.defit, order=c(1, 0, 1), include.mean=F, method = "CSS-ML")
summary(goog.fit)
goog.pred = predict(goog.fit, n.ahead=100)


plot(c(goog.defit[2600:n], goog.pred$pred), type="l", lty=2)
lines(goog.defit[2600:n], type="l")
