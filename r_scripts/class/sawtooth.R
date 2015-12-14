sawtooth <- function(ncycles = 5, amp = 10, plot.it = as.logical(T)) {
#----------------------------------------------------------------------
#  R function to generate and plot a data set consisting of ncycles
#  cycles of the sawtooth function with amplitude amp.
#
#  The function sawtooth has three optional input arguments:
#   ncycles = an integer scalar indicating the number of cycles 
#             to generate. The default is ncycles = 5.
#   amp = an integer scalar containing the amplitude of the cycles.
#   plot.it = a logical indicating whether or not to generate the 
#             plot.  The default is plot.it = TRUE.
#
#  The function sawtooth function returns a list containing two objects
#    t = an integer vector containing the abscissa of the points
#    x = an real vector containing the ordinates of the points from 
#        the sawtooth function having absciss t.
#  If plot.it = TRUE, then the sawtooth function will also create a
#  graph of the sawtooth function having ncycles cycles and amplitude
#  as specified in amp.
#----------------------------------------------------------------------
  n <- 4*amp*ncycles
  t <- 1:n
  x <- rep(0, n)
  a2 <- 2*amp
#
  x[1:amp]       <- t[1:amp]
  x[(n-amp+1):n] <- -amp + t[1:amp]
#
  for(cycle in 1:(2*ncycles-1)) {
      cc <- 2*cycle - 1
      dd <- 2*cycle + 1
      x[(cc*amp + 1):(dd*amp)] <- ((-1)^(cycle+1))*amp + ((-1)^cycle)*t[1:a2]
  }
#
  if(plot.it)
    plot(x, type = "l", xlab = "t", ylab = expression(x[t]),
         main = "Sawtooth Function")

  return(list(t = t, x = x))

}