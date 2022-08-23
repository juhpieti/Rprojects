library(laGP) # for rmvn function to sample multivariate normal

source("model_metrics.R")

### assume that we have fitted our model with test points x and have 2000 test points y
### and that we also have the predictive distribution y|x where x is the training data
### we will then simulate N predictions from that distribution for every y_i in the test set
### let's call y*_ij the j:th simulated prediction

### we want to end up with vector u so that
### u_i = (1/n)*sum_j(I(y*_ij < y_i))

### so we are taking the proportions of predictions y*_ij having lower value than the observed y_i
### that vector u should have an Uniform(0,1) distribution

### y = vector of realized observations (not used to build the model)
### pred_means = vector of corresponding means from predictive distribution y|x
### pred_covmat = predictive covariance matrix
### pred_sds = vector of corresponding standard deviations from predictive distribution
### m = number of posterior predictive samples

PIT_test <- function(y, pred_means, pred_covmat, m = 250, return_predictions = TRUE) {
  
  n <- length(y)
  u <- rep(0, n)
  
  pred_sample <- rmvn(m, pred_means, pred_covmat) # samples vectors from predictive mvnm
  
  for (i in 1:n) {
    u[i] <- mean(pred_sample[,i] < y[i])
  }
  
  hist(u, breaks = 9)
  if (return_predictions) {
    return(pred_sample)
  }
}

### I'm doing another function too, that doesn't take covariance matrix but just the diagonal of it = variances
### just because at least predict.hetGP() gives you only the predictive variances 

PIT_test_sd <- function(y, pred_means, pred_sds, m = 250, return_predictions = TRUE) { 
  
  n <- length(y)
  u <- rep(0, n)
  
  pred_sample <- matrix(0, m, n)
  
  for (i in 1:n) {
    pred_sample[,i] <- rnorm(m, pred_means[i], pred_sds[i]) # samples values from predictive marginal univariate normals
    u[i] <- mean(pred_sample[,i] < y[i])
  }
  
  hist(u, breaks = 9)
  
  if (return_predictions) {
    return(pred_sample)
  }
}

##### experimenting #####

### let's try with i.i.d. sample from univariate N(0,1)

n <- 2000

y <- rnorm(n, 0, 1)
hist(y)

sample <- PIT_test_sd(y, rep(0,n),rep(1,n), m = 100) # looks like uniformly distributed, fine!

### let's test what happens if we generate some data from known multivariate normal distribution
### and then sample from that known mvnm (more related to our task with GP's?)

n <- 2000

x <- seq(0,100, length.out = n)
X <- distance(x) # function from package "laGP"
X <- exp(-X*0.25)
y <- rmvn(1, rep(0, length(x)), X)

plot(x,y, type = "l")

sample1 <- PIT_test_sd(y, rep(0, length(x)), sqrt(diag(X)), m = 100) # not so promising as the previous
sample2 <- PIT_test(y, rep(0, length(x)), X, m = 100) 


### trying to illustrate how the u vector is constructed (samples from marginal univariate normals)

plot(NULL,NULL, ylim = c(-5,5), xlim = c(0,100))

for (i in 1:length(x)) {
  for (j in 1:nrow(sample1)) {
    points(x = x[i], y = sample1[j,i], col = 'red', cex = .3)
  }
}

points(x,y, pch = 16, cex = .5, type = "l")

### with the second approach too (samples from mvnm)

plot(NULL,NULL, ylim = c(-5,5), xlim = c(0,100))

for (i in 1:length(x)) {
  for (j in 1:nrow(sample2)) {
    points(x = x[i], y = sample2[j,i], col = 'red', cex = .3)
  }
}

points(x,y, pch = 16, cex = .5, type = "l")


##### from the following I got some pictures of how the histograms looked like with hetGP & mgcv #####

### test with hetGP model

test_data <- shuffle_data(des_mat_10_10K)

pred <- fit_model(n_knots = 1000, n_test_points = 2000, pred_type = "original", package = "hetGP",
                      nugget_known = FALSE, known_data = test_data, verb = 1)

PIT_test_sd(pred$observed, pred$pred, pred$pred_se, m = 250, return_predictions = FALSE)

### test with mgcv model

pred1 <- fit_model(n_knots = 1000, n_test_points = 2000, pred_type = "original", package = "mgcv", n_interactions = 5,
                   known_data = test_data, verb = 1)

PIT_test_sd(pred1$observed, pred1$pred, pred1$pred_se, m = 250, return_predictions = FALSE)

### test with laGP model

pred2 <- fit_model(n_knots = 1000, n_test_points = 2000, pred_type = "normal", package = "laGP",
                       nugget_known = FALSE, known_data = test_data, verb = 1)

PIT_test_sd(pred2$observed, pred2$pred, pred2$pred_se, m = 250, return_predictions = FALSE)

# PIT_test(y, pred2$mean, pred2$Sigma, m = 100, return_predictions = FALSE) # fit_model() needs to be modified to return cov matrix for this
