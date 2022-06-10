n <- 100
X <- matrix(seq(0, 10, length=n), ncol=1)

library(plgp)
D <- distance(X) #eucledian distances ||x-x'||^2 

eps <- sqrt(.Machine$double.eps)
Sigma <- exp(-D) + diag(eps, n) #some jitter on diagonal to prevent zeros?

library(mvtnorm)
Y <- rmvnorm(1, sigma=Sigma)

plot(X, Y, type="l")

Y <- rmvnorm(5, sigma=Sigma)
matplot(X, t(Y), type="l", ylab="Y")
abline(h = c(-2,2), col = 'red', type = "l", lty = 2, lwd = 3)

n <- 8
X <- matrix(seq(0,2*pi,length=n), ncol=1) #training points
y <- sin(X)
D <- distance(X)
Sigma <- exp(-D) + diag(eps, ncol(D))

XX <- matrix(seq(-0.5, 2*pi+0.5, length=100), ncol=1) #test points
DXX <- distance(XX)
SXX <- exp(-DXX) + diag(eps, ncol(DXX))

DX <- distance(XX, X)
SX <- exp(-DX)

Si <- solve(Sigma)
mup <- SX %*% Si %*% y
Sigmap <- SXX - SX %*% Si %*% t(SX)

YY <- rmvnorm(100, mup, Sigmap) ## sampling from posterior with 8 obsevrations y
q1 <- mup + qnorm(0.05, 0, sqrt(diag(Sigmap)))
q2 <- mup + qnorm(0.95, 0, sqrt(diag(Sigmap)))


### this is interesting one
matplot(XX, t(YY), type="l", col="gray", lty=1, xlab="x", ylab="y")
points(X, y, pch=20, cex=2)
lines(XX, sin(XX), col="blue")
lines(XX, mup, lwd=2)
lines(XX, q1, lwd=2, lty=2, col=2)
lines(XX, q2, lwd=2, lty=2, col=2)


