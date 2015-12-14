source('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab.R')
dyn.load('C:/Dropbox/Baylor/TimeSeries/R/ts_hw/ts_hw/timeslab2015_64.dll')

perdave(n=10, alpha=c(0.05, 0.09), nsamps = 20, Q = 10, rvar = 1, seed = 0)
arsp(c(0.05, 0.09), 1, 10)


### ADV COMP ORG EXAM #2
log4 = function(x){return(log2(x)/log2(4))}
entropy = function(x){return(x*log2(1/x))}
entropy4 = function(x){return(x*log4(1/x))}
sysentropy = function(x){
	sum = 0;
	for (i in x) { sum = sum + entropy(i)}
	return(sum)
}
sysentropy4 = function(x){
	sum = 0;
	for (i in x) { sum = sum + entropy4(i)}
	return(sum)
}