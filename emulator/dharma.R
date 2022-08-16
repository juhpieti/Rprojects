library(DHARMa)
source("model_metrics.R")

### y = observed values
### pred_means = predictions = posterior means
### pred_sds = prediction standard deviation (dharma_plots_cov() will take the whole predictive covariance matrix)
### m = number of simulations

### plots you DHARMa-plots and returns DHARMa object
### simulations/samples to create DHARMa object are done by sampling from univariate normal distributions with
### given means and sd's 

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

### similar to dharma_plots() but instead to create DHARMa-object,
### it samples from multivariate normal using covariance matrix given
dharma_plots_cov <- function(y, pred_means, cov_matrix, m = 250) {
  
  n <- length(y)
  fitted <- pred_means
  pred_sample <- rmvn(m, fitted,cov_matrix) # samples from predictive distribution (multivariate normal)
  
  dharma <- createDHARMa(t(pred_sample), y, fitted)
  plot(dharma)
  return(dharma)
}


##### tests #####

### laGP, hetGP, mgcv

testing_data <- shuffle_data(des_mat_10_10K) # set the data (training set, test set to experiment with)

lagpmod <- fit_model(n_knots = 1000, package = "laGP", pred_type = "quantile", nugget_known = FALSE, known_data = testing_data, verb = 1)
dharma_lagp <- dharma_plots(lagpmod$observed, lagpmod$pred, lagpmod$pred_se)

# lagpmod_cov <- fit_model(n_knots = 1000, n_test_points = 2000, package = "laGP", pred_type = "quantile", nugget_known = FALSE, known_data = testing_data)
# dharma_plots_cov(real_obs, lagpmod_cov$mean, lagpmod_cov$Sigma) # requires tuning of fit_model() to return cov_matrix

hetmod <- fit_model(n_knots = 1000, package = "hetGP",pred_type = "normal", nugget_known = FALSE, known_data = testing_data, verb = 1)
dharma_het <- dharma_plots(hetmod$observed, hetmod$pred, hetmod$pred_se, m = 250)

testDispersion(dharma_het)
testUniformity(dharma_het)
testQuantiles(dharma_het)
plotResiduals(dharma_het)

mgcvmod <- fit_model(n_knots = 1000, package = "mgcv", gam_interact = 5, known_data = testing_data, verb = 1)
dharma_mgcv <- dharma_plots(mgcvmod$observed, mgcvmod$pred, mgcvmod$pred_se, m = 250)

# ### just for curiosity with known, univariate normal, does it look perfect?
# x <- rep(1:20,100)
# y <- rnorm(2000, 3*x, 1) + 5
# 
# dharma_plots(y, 3*x + 5, rep(1, 10000), m = 250)
# 
# ### yes it does!


# ### mgcv is included in DHARMa "supported packages" so you actually just need the model and DHARMa does the rest
# 
# mgcvmodel <- fit_model(n_knots = 1000, n_test_points = 2000, package = "mgcv", n_interactions = 5, known_data = testing_data)
# simulationOutput <- simulateResiduals(fittedModel = mgcvmodel, plot = T)
# plotResiduals(simulationOutput)