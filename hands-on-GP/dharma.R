library(DHARMa)
source("model_metrics.R")

### y = observed values
### pred_means = predictions = posterior means
### pred_sds = prediction standard deviation (dharma_plots_cov() will take the whole predictive covariance matrix)
### m = number of simulations

dharma_plots <- function(y, pred_means, pred_sds, m = 250) { 
  
  n <- length(y)
  fitted <- pred_means
  pred_sample <- matrix(0, n, m)
  
  for (i in 1:n) {
    pred_sample[i,] <- rnorm(m, pred_means[i], pred_sds[i]) # samples values from predictive marginal univariate normals
  }
  
  dharma <- createDHARMa(pred_sample, y, fitted)
  plot(dharma)
  return(dharma)
}

dharma_plots_cov <- function(y, pred_means, cov_matrix, m = 250) {
  
  n <- length(y)
  fitted <- pred_means
  pred_sample <- rmvn(m, fitted,cov_matrix) # samples from predictive distribution (multivariate normal)
  
  dharma <- createDHARMa(t(pred_sample), y, fitted)
  plot(dharma)
  return(dharma)
}

##### tests #####

testing_data <- shuffle_data(des_mat_10_10K)
real_obs <- testing_data$likelihood[1:2000]

lagpmod <- fit_model(n_knots = 1000, n_test_points = 2000, package = "laGP", pred_type = "quantile", nugget_known = FALSE, known_data = testing_data, verb = 1)
dharma_lagp <- dharma_plots(real_obs, lagpmod$pred, lagpmod$pred_se)

# lagpmod_cov <- fit_model(n_knots = 1000, n_test_points = 2000, package = "laGP", pred_type = "quantile", nugget_known = FALSE, known_data = testing_data)
# dharma_plots_cov(real_obs, lagpmod_cov$mean, lagpmod_cov$Sigma) # requires tuning of fit_model() to return cov_matrix

hetmod <- fit_model(n_knots = 1000, n_test_points = 2000, package = "hetGP",pred_type = "normal", nugget_known = FALSE, known_data = testing_data, verb = 1)
dharma_obj <- dharma_plots(real_obs, hetmod$pred, hetmod$pred_se, m = 250)

testDispersion(dharma_obj)
testUniformity(dharma_obj)
testQuantiles(dharma_obj)
plotResiduals(dharma_obj)

mgcvmod <- fit_model(n_knots = 1000, n_test_points = 2000, package = "mgcv", n_interactions = 5, known_data = testing_data, verb = 1)
dharma_plots(real_obs, mgcvmod$pred, mgcvmod$pred_se, m = 250)

# 
# y <- rnorm(1000, 0, 1)
# 
# dharma_plots(y, rep(0, 1000), rep(1, 1000), m = 100)

mgcvmodel <- fit_model(n_knots = 1000, n_test_points = 2000, package = "mgcv", n_interactions = 5, known_data = testing_data)
simulationOutput <- simulateResiduals(fittedModel = mgcvmodel, plot = T)
plotResiduals(simulationOutput)

hetmod
