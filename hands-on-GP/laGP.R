library(laGP)
library(tidyverse)
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
load("metrics_matrix_new.Rdata")
load("priorfcn2Juho.Rdata")
source("model_metrics.R")

deleteGPsep(laGPmodel)
deleteGPsep(laGPmodel_01)

shuffled_df <- model_metrics(1, 1, 100, 150, "mlegp")

mod <- model_metrics(n_knots = 600, package = "laGP")

for (i in 1:10) {
  model_metrics(site_no = 1, response = 1, n_knots = 600, n_test_points = 150, package = "laGP", predictor_type = "quantile")
}

model_metrics(site_no = 1, response = 1, n_knots = 600, n_test_points = 150, package = "laGP", predictor_type = "quantile")

df <- as.data.frame(SS.stack[[1]][[1]])
colnames(df)[6] <- "likelihood"
shuffled_df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
shuffled_df_01 <- shuffled_df

shuffled_df_01$som_respiration_rate <- sapply(shuffled_df_01$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
shuffled_df_01$rdConst <- sapply(shuffled_df_01$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
shuffled_df_01$psnTOpt <- sapply(shuffled_df_01$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
shuffled_df_01$wueConst <- sapply(shuffled_df_01$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
shuffled_df_01$leafGrowth <- sapply(shuffled_df_01$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})

training_set <- shuffled_df[(151:250), ]
training_set <- shuffled_df[1:100, ] ######
training_set_01 <- shuffled_df_01[(151:250), ]


test_set <- shuffled_df[1:150, ]
test_set <- shuffled_df[1:750, ] #####
test_set_01 <- shuffled_df_01[1:150, ]

X = training_set[, -ncol(df), drop = FALSE]
X_01 = training_set_01[, -ncol(df), drop = FALSE]

Z = training_set[, ncol(df), drop = TRUE]
Z_01 = training_set_01[, ncol(df), drop = TRUE]


d <- darg(list(mle = TRUE, min = 0), X)
d_01 <- darg(list(mle = TRUE, min = 0), X_01)



#d <- darg(list(mle=FALSE), X)
g <- garg(list(mle=FALSE), Z) #mle=FALSE = no estimation for nugget term, just set to 0
g_01 <- garg(list(mle=FALSE), Z_01) 


laGPmodel <- newGPsep(X = X, Z = Z, d = rep(d$start, 5), g = 0, dK = TRUE)
laGPmodel_01 <- newGPsep(X = X_01, Z = Z_01, d = rep(d_01$start, 5), g = 0, dK = TRUE)

#laGPmodel_01 <- newGPsep(X = X_01, Z = Z_01, d = c(1/1.201387381, 1/0.002795255, 1/5.813077585, 1/32.433430802, 1/4.678862480), g = 0, dK = TRUE)




#laGPmodel <- newGPsep(X = X, Z = Z,
#                      d = c(4.550509e+00, 7.776160e-04, 1.167631e-03, 6.664508e-03, 1.362268e-05), g = 0, dK = FALSE)


#mleGPsep(laGPmodel, param="d", tmin=c(d$min, 0), tmax=c(d$max, 0), verb = 1)
mleGPsep(laGPmodel, param="d", tmin=d$min, tmax=d$max, verb = 1)

mleGPsep(laGPmodel_01, param="d", tmin = d_01$min, tmax = d_01$max, verb = 1)

#mleGPsep(laGPmodel_01, param="d", tmin=c(d_01$min, 0), tmax=c(d_01$max, 0), verb = 1)
#mleGPsep(laGPmodel_01, param="d", tmin = d_01$min, tmax = d_01$max*1e05, ab = d_01$ab, verb = 1)
#mleGPsep(laGPmodel_01, param="d", tmin = rep(d_01$min, 5), tmax = rep(d_01$max, 5), verb = 1)
#mleGPsep(laGPmodel_01, param="d", tmin = .Machine$double.eps, tmax = d_01$max, verb = 1)

#model_metrics(site_no = 1, response = 1, n_knots = 100, n_test_points = 150, predictor_type = "quantile", package = "mlegp", data = shuffled_df) # mlegp!!!


laGPpred <- predGPsep(laGPmodel, test_set[, -ncol(df), drop = FALSE], lite = TRUE, nonug = FALSE)
laGPpred_01 <- predGPsep(laGPmodel_01, test_set_01[, -ncol(df), drop = FALSE], lite = TRUE, nonug = FALSE)


test_set_mod <- test_set
test_set_mod$pred = laGPpred$mean
test_set_mod$se = sqrt(laGPpred$s2)

sqrt(mean((test_set_mod$likelihood - test_set_mod$pred)^2)) #1135
mean(test_set_mod$se) #1194

test_set_mod_01 <- test_set_01
test_set_mod_01$pred = laGPpred_01$mean
test_set_mod_01$se = sqrt(laGPpred_01$s2)

sqrt(mean((test_set_mod_01$likelihood - test_set_mod_01$pred)^2)) #1135
mean(test_set_mod_01$se) #1194

deleteGPsep(laGPmodel)
deleteGPsep(laGPmodel_01)

plot(test_set_mod$likelihood, test_set_mod$pred, ylim = c(6000, 13000))
abline(0, 1)


##### loops for running laGP #####

### on Thursday run a couple more for site 1, largest no_knots
for (n in c(100, 200, 300, 400, 500, 600)) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = 1, response = resp, n_knots = n, n_test_points = 150, package = "laGP", predictor_type = "normal")
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(1, resp, n, 150, as.numeric(output[[1]]), as.numeric(output[[2]]), as.numeric(output[[3]]), as.numeric(output[[4]]), "laGP", "normal")
  }
}

model_metrics(package = "laGP")

###for different sites with n_knots = 100
for (site in 1:12) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = site, response = resp, n_knots = 600, n_test_points = 150, package = "laGP", predictor_type = "normal")
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(site, resp, 100, 150, as.numeric(output[[1]]), as.numeric(output[[2]]), as.numeric(output[[3]]), as.numeric(output[[4]]), "laGP", "normal")
    #print(model_metrics(site_no = site, response = resp, n_knots = 600, n_test_points = 150, package = "laGP", predictor_type = "normal"))
  }
} 

for (n in c(100, 200, 300, 400, 500, 600)) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = 1, response = resp, n_knots = n, n_test_points = 150, package = "laGP", predictor_type = "normal")
    print(output)
    #metrics_matrix_new[nrow(metrics_matrix_new)+1, ] <- c(1, resp, n, 150, output[[1]], output[[2]], output[[3]], output[[4]], "laGP")
  }
}

### testing if there's difference between using g = 0, g = .Machine$double.eps
rmse_normal <- rep(0, 100)
rmse_quantile <- rep(0,100)
mstde_normal <- rep(0,100)
mstde_quantile <- rep(0,100)
time_normal <- rep(0,100)
time_quantile <- rep(0,100)

for (i in 1:100) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_normal <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
                                  package = "laGP", predictor_type = "normal", data = df)
  metrics_quantile <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
                                  package = "laGP", predictor_type = "quantile", data = df)
  rmse_normal[i] <- metrics_normal[[3]]
  rmse_quantile[i] <- metrics_quantile[[3]]
  mstde_normal[i] <- metrics_normal[[4]]
  mstde_quantile[i] <- metrics_quantile[[4]]
  time_normal[i] <- metrics_normal[[1]] + metrics_normal[[2]]
  time_quantile[i] <- metrics_quantile[[1]] + metrics_quantile[[2]]
}

df_metrics <- data.frame(rmse_normal = rmse_normal, rmse_quantile = rmse_quantile, mstde_normal = mstde_normal,
                         mstde_quantile = mstde_quantile, time_normal = time_normal, time_quantile = time_quantile)

save(metrics_matrix_new, file = "metrics_matrix_new.Rdata")
save(metrics_matrix, file = "metrics_matrix.Rdata")

plot21 <- metrics_matrix_new %>%
  filter(site == 1, package == "laGP") %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color=knots, shape=response)) +
  geom_point(size = 3) + labs(x = "time", title = "Effects of increasing number of knots using site No.1", subtitle = "package: laGP")

plot22 <- metrics_matrix_new %>%
  filter(site == 1, package == "laGP") %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color=knots, shape=response)) +
  geom_point(size = 3) + labs(x = "time")

plot11 <- metrics_matrix_new %>%
  filter(response == 1, site == 1) %>%
  ggplot(mapping = aes(x = time_model + time_pred, y = RMSqE, color = knots, shape = package)) +
  geom_point(size = 3) + labs(x = "time", title = "Differences between mlegp and laGP", subtitle = "site No. = 1, output = 1")

plot12 <- metrics_matrix_new %>%
  filter(response == 1, site == 1) %>%
  ggplot(mapping = aes(x = time_model + time_pred, y = MStdE, color = knots, shape = package)) +
  geom_point(size = 3) + labs(x = "time", title = "Differences between mlegp and laGP", subtitle = "site No. = 1, output = 1")
  


(plot1 <- ggpubr::ggarrange(plot11, plot12, nrow = 2, common.legend = TRUE, legend = "right"))
(plot2 <- ggpubr::ggarrange(plot21, plot22, nrow = 2, common.legend = TRUE, legend = "right"))

metrics_matrix_new %>%
  filter(package == "laGP") %>%
  ggplot(mappin = aes(x = time_model + time_pred, y = RMSqE, color = knots, shape = response)) +
  geom_point()

metrics_matrix_new %>%
  filter(package == "laGP") %>%
  ggplot(mappin = aes(x = time_model + time_pred, y = MStdE, color = knots, shape = response)) +
  geom_point()
  

tic <- proc.time()
GPmodel <- mlegp(X = training_set[, -ncol(df), drop = FALSE], Z = training_set[, ncol(df), drop = FALSE], 
                 nugget = 0, nugget.known = 1, verbose = 0)
toc <- proc.time()
time_model = as.numeric((toc - tic)[3])

tic <- proc.time()
GP_pred = predict.gp(GPmodel, test_set[ ,-ncol(df), drop = FALSE], se.fit = TRUE)

toc <- proc.time()
time_pred <- as.numeric((toc - tic)[3])


test_set <- cbind(test_set, GP_pred$fit, GP_pred$se.fit)
colnames(test_set)[7:ncol(test_set)] <- c("pred", "se")

rmse <- sqrt(mean((test_set$likelihood - test_set$pred)^2, na.rm = TRUE)) 
mean_std_err <- mean(test_set$se, na.rm = TRUE)

metrics <- c("model time" = time_model, "prediction time" = time_pred, "rooted mean square error" = rmse,
             "mean standard error" = mean_std_err)

plot(GPmodel)
return(metrics)



##### experimenting with "normal" predictors vs quantile predictors #####

rmse_normal <- rep(0, 10)
rmse_quantile <- rep(0,10)
mstde_normal <- rep(0,10)
mstde_quantile <- rep(0,10)
time_normal <- rep(0,10)
time_quantile <- rep(0,10)

for (i in 1:10) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_normal <- model_metrics(site_no = 1, response = 1, n_knots = 150, n_test_points = 150, package = "laGP", predictor_type = "normal", data = df)
  metrics_quantile <- model_metrics(site_no = 1, response = 1, n_knots = 150, n_test_points = 150, package = "laGP", predictor_type = "quantile", data = df)
  rmse_normal[i] <- metrics_normal[[3]]
  rmse_quantile[i] <- metrics_quantile[[3]]
  mstde_normal[i] <- metrics_normal[[4]]
  mstde_quantile[i] <- metrics_quantile[[4]]
  time_normal[i] <- metrics_normal[[1]] + metrics_normal[[2]]
  time_quantile[i] <- metrics_quantile[[1]] + metrics_quantile[[2]]
}

df_metrics <- data.frame(rmse_normal = rmse_normal, rmse_quantile = rmse_quantile, mstde_normal = mstde_normal,
                         mstde_quantile = mstde_quantile, time_normal = time_normal, time_quantile = time_quantile)
  
model_metrics(site_no = 1, response = 1, n_knots = 150, n_test_points = 150, package = "laGP", predictor_type = "quantile")

##### EXAMPLES #####

X <- matrix(seq(0,2*pi,length=7), ncol=1)
Z <- sin(X)

## new GP fit
gpi <- newGP(X, Z, 2, 0.1, dK = TRUE)
#gpi <- newGP(X, Z, 1, 0.000001)
#gpi <- newGP(X, Z, 2, 0.01)

## make predictions
XX <- matrix(seq(-1,2*pi+1, length=499), ncol=ncol(X))

d <- darg(NULL, X)
g <- garg(list(mle=TRUE), Z)

mleGP(gpi, param="both", tmin=c(d$min, g$min), tmax=c(d$max, g$max))
print(jmleGP(gpi, drange=c(d$min, d$max), grange=c(g$min, g$max)))


p <- predGP(gpi, XX, nonug = TRUE)

deleteGP(gpi)

## sample from the predictive distribution
library(mvtnorm)
N <- 100
ZZ <- rmvt(N, p$Sigma, p$df)
ZZ <- ZZ + t(matrix(rep(p$mean, N), ncol=N))
matplot(XX, t(ZZ), col="gray", lwd=0.5, lty=1, type="l", 
        xlab="x", ylab="f-hat(x)", bty="n")
points(X, Z, pch=19, cex=2)

## update with four more points
X2 <- matrix(c(pi/2, 3*pi/2, -0.5, 2*pi+0.5), ncol=1)
Z2 <- sin(X2)
updateGP(gpi, X2, Z2)

## make a new set of predictions
p2 <- predGP(gpi, XX)
ZZ <- rmvt(N, p2$Sigma, p2$df) 
ZZ <- ZZ + t(matrix(rep(p2$mean, N), ncol=N))
matplot(XX, t(ZZ), col="gray", lwd=0.5, lty=1, type="l", 
        xlab="x", ylab="f-hat(x)", bty="n")
points(X, Z, pch=19, cex=2)
points(X2, Z2, pch=19, cex=2, col=2)

## clean up
deleteGP(gpi)

## a simple example with estimated nugget
library(MASS)

## motorcycle data and predictive locations
X <- matrix(mcycle[,1], ncol=1)
Z <- mcycle[,2]

## get sensible ranges
d <- darg(NULL, X)
g <- garg(list(mle=TRUE), Z)

## initialize the model
gpi <- newGP(X, Z, d=d$start, g=g$start, dK=TRUE)

## separate marginal inference (not necessary - for illustration only)
print(mleGP(gpi, "d", d$min, d$max))
print(mleGP(gpi, "g", g$min, g$max))

## joint inference (could skip straight to here)
print(jmleGP(gpi, drange=c(d$min, d$max), grange=c(g$min, g$max)))

## plot the resulting predictive surface
N <- 100
XX <- matrix(seq(min(X), max(X), length=N), ncol=1)
p <- predGP(gpi, XX, lite=TRUE)
plot(X, Z, main="stationary GP fit to motorcycle data")
lines(XX, p$mean, lwd=2)
lines(XX, p$mean+1.96*sqrt(p$s2*p$df/(p$df-2)), col=2, lty=2)
lines(XX, p$mean-1.96*sqrt(p$s2*p$df/(p$df-2)), col=2, lty=2)

## clean up
deleteGP(gpi)

## 
## with a separable correlation function 
##

## 2D Example: GoldPrice Function, mimicking GP_fit package
f <- function(x) 
{
  x1 <- 4*x[,1] - 2
  x2 <- 4*x[,2] - 2;
  t1 <- 1 + (x1 + x2 + 1)^2*(19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2);
  t2 <- 30 + (2*x1 -3*x2)^2*(18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2);
  y <- t1*t2;
  return(y)
}

## build design
library(tgp)
n <- 100 ## change to 100 or 1000 for more interesting demo
B <- rbind(c(0,1), c(0,1))
X <- dopt.gp(n, Xcand=lhs(10*n, B))$XX
## this differs from GP_fit in that we use the log response
Y <- log(f(X))

## get sensible ranges
d <- darg(NULL, X)
g <- garg(list(mle=TRUE), Y)

## build GP and jointly optimize via profile mehtods
gpisep <- newGPsep(X, Y, d=rep(d$start, 2), g=g$start, dK=TRUE)
jmleGPsep(gpisep, drange=c(d$min, d$max), grange=c(g$min, g$max))

## clean up
deleteGPsep(gpisep)

## alternatively, we can optimize via a combined gradient
gpisep <- newGPsep(X, Y, d=rep(d$start, 2), g=g$start, dK=TRUE)
mleGPsep(gpisep, param="both", tmin=c(d$min, g$min), tmax=c(d$max, g$max))
deleteGPsep(gpisep)


## a "computer experiment" -- a much smaller version than the one shown
## in the aGP docs

## Simple 2-d test function used in Gramacy & Apley (2015);
## thanks to Lee, Gramacy, Taddy, and others who have used it before
f2d <- function(x, y=NULL)
{
  if(is.null(y)) {
    if(!is.matrix(x) && !is.data.frame(x)) x <- matrix(x, ncol=2)
    y <- x[,2]; x <- x[,1]
  }
  g <- function(z)
    return(exp(-(z-1)^2) + exp(-0.8*(z+1)^2) - 0.05*sin(8*(z+0.1)))
  z <- -g(x)*g(y)
}

## design with N=441
x <- seq(-2, 2, length=11)
X <- expand.grid(x, x)
Z <- f2d(X)

## fit a GP
gpi <- newGP(X, Z, d=0.35, g=1/1000)

## predictive grid with NN=400
xx <- seq(-1.9, 1.9, length=20)
XX <- expand.grid(xx, xx)
ZZ <- f2d(XX)

## predict
p <- predGP(gpi, XX, lite=TRUE)
## RMSE: compare to similar experiment in aGP docs
sqrt(mean((p$mean - ZZ)^2))

## visualize the result
par(mfrow=c(1,2))
image(xx, xx, matrix(p$mean, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="predictive mean")
image(xx, xx, matrix(p$mean-ZZ, nrow=length(xx)), col=heat.colors(128),
      xlab="x1", ylab="x2", main="bas")

## clean up
deleteGP(gpi)

## see the newGP and mleGP docs for examples using lite = FALSE for
## sampling from the joint predictive distribution
