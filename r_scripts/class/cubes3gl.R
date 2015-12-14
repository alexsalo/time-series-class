library('rgl')
install.packages('magrittr')
library('magrittr')

cube <- function(x=0,y=0,z=0, bordered=TRUE, 
                 filled = TRUE, lwd=2, scale=1,
                 fillcol = gray(.95),
                 bordercol ='black', ...) {
  
  mycube <- cube3d()
  
  # Reduce size to unit
  mycube$vb[4,] <- mycube$vb[4,]/scale*2
  
  for (i in 1:length(x)) {
    # Add cube border
    if (bordered) {
      bcube <- mycube
      bcube$material$lwd <- lwd
      bcube$material$front <- 'line'
      bcube$material$back <- 'line'
      bcube %>% translate3d(x[i], y[i], z[i]) %>% shade3d
    }
    # Add cube fill
    if (filled) {
      fcube <- mycube
      fcube$vb[4,] <- fcube$vb[4,]*1.01
      fcube$material$col <- fillcol
      fcube %>% translate3d(x[i], y[i], z[i]) %>% shade3d
    }
  }
}

clear3d()
cube(0,0,0)
cube(1,0,0, filled=F)
cube(-1,0,0, bordered=F)
movie3d(spin3d(axis=c(0,0,1), rpm=20), duration=2.95) 


# I mapped R using an excel spreadsheet which 
# translated Xs into 2D locations points
clear3d()
y <- c(1,1,1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,
       5,5,6,6,6,6,7,7,7,7,7,8,8,8,8,9,9,9,9,10,
       10,10,10,11,11,11,11,11,12,12,12,12,12)

x <- c(8,7,6,3,2,1,7,6,3,2,6,5,3,2,6,5,4,3,2,5,4,
       3,2,5,4,3,2,6,5,4,3,2,7,6,3,2,7,6,3,2,7,6,
       3,2,6,5,4,3,2,5,4,3,2,1)

z <- cummax(y)*.5

length(x)==length(y)
cube(x,y,z)
movie3d(spin3d(axis=c(0,0,1), rpm=20), duration=2.95) 


# Let's see how sin and cos can work together
z <- seq(0,6,.1)
x <- sin(pi*z)*z
y <- cos(pi*z)*z

clear3d()
cube(x,y,z*2, scale=.75)
movie3d(spin3d(axis=c(0,0,1), rpm=20), duration=2.95) 