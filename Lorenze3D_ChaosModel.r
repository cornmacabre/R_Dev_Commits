# Some parameter values
w <- 28.0
s <- 10.0
b <- 2.66666666
vars <- 3 # The number of variables in the system.
# Some initial conditions
x <- c(5,5,5)
h <- 0.01 # The step size for the Runge-Kutta
timeserieslength <- 2024 # The number of iterations we will use for the Runge-Kutta
# Initializing some vectors.
deriv <- 1:vars;m <- x;xNEW <- x;
xk1 <- 1:vars;xk2 <- 1:vars;xk3 <- 1:vars;xk4 <- 1:vars;

# We define the Lorenz equations as an R function.
eqs <- function(n){
deriv[1] <- s*(n[2] - n[1])
deriv[2] <- (w*n[1]) - n[2] - (n[1]*n[3])
deriv[3] <- (n[1]*n[2]) - (b*n[3])
return(deriv)
}

# Now we enter the fourth-order Runge-Kutta.
xe <- NULL
for (i in 1:timeserieslength) {
m <- x
derivatives <- eqs(m)
xk1 <- derivatives
m <- x+(.5*h*xk1)
derivatives <- eqs(m)
xk2 <- derivatives
m <- x+(.5*h*xk2)
derivatives <- eqs(m)
xk3 <- derivatives
m <- x + h*xk3
derivatives <- eqs(m)
xk4 <- derivatives
xNEW <- x+((h/6)*(xk1 + (2*xk2) + (2*xk3) + xk4));
# Now we need to collect our computed values.
xe <- rbind(xe,x)

# Now we assign the new variable values to the current values to do the next RK4 iteration.
x <- xNEW
}
# We are done with the Runge-Kutta.

cxe <- cbind(xe,c(1:nrow(xe))) # Adding a counting variable to the data matrix
# cxe

# Now we plot the Lorenz attractor.
# First we get the endpoint coordinates for the axes.
x1max <- max(cxe[,1]);x1min <- min(cxe[,1])
x2max <- max(cxe[,2]);x2min <- min(cxe[,2])
x3max <- max(cxe[,3]);x3min <- min(cxe[,3])
origin <- cbind(x1min,x2min,x3min,timeserieslength+2)
x1next <- cbind(x1max,x2min,x3min,timeserieslength+3)
x2next <- cbind(x1min,x2max,x3min,timeserieslength+4)
x3next <- cbind(x1min,x2min,x3max,timeserieslength+5)
# We are done getting the endpoint coordinates for the axes.
final <- cbind(NA,NA,NA,timeserieslength+1) #To break the plotting of axis lines from the other data
# Now we bind the axes coordinates to the original data matrix
axes <- rbind(origin,rbind(x1next,rbind(origin,rbind(x2next,(rbind(origin,rbind(x3next,final)))))))
trajects <- rbind(axes,cxe)
# trajects
# Translating a 3-D projection onto a 2-D plot
angle <- 2*pi/3
xvalue <- trajects[,1] - ((cos(angle))*trajects[,3])
yvalue <- trajects[,2] - ((sin(angle))*trajects[,3])
# The final plot statement for display
plot(xvalue, yvalue, axes=FALSE, xlab="", ylab="", type="l", lwd=2)

# A pdf version of the plot saved to the R workspace
savePlot(filename = "Lorenz",
type = "pdf")

# Now a 3-D display of the attractor that is rotatable
library(rgl) 
plot3d(cxe,lwd=4,type="l",col=rainbow(1000))