library(DHARMa)

### bunch of functions used to simplify the major script model_metrics.R

### takes in design_matrix, shuffles and RETURNS it
### e.g. SS.stack has 12 lists (different sites) of 2 data frames (response 1 and 2) so you have to specify site and response to use
### can be used to create a common data to use with different settings / packages to compare performance
shuffle_data <- function(design_matrix = SS.stack, site_no = 1, response = 1) {
  ## check whether your data is just 1 resp 1 site or [[]][[]] type (e.g. SS.stack)?
  if (length(design_matrix[[1]][[1]]) > 1) { #SS.stack with 12 sites, 2 responses
    df <- as.data.frame(design_matrix[[site_no]][[response]])
  } else {
    df <- as.data.frame(design_matrix)
  }
  colnames(df)[ncol(df)] <- "likelihood"
  df_shuffled <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ] # shuffles the data
  
  return(df_shuffled)
}


### RETURNS shuffled and prepared (modifies parameter space) design_matrix given to use for modeling
### if known_data given (e.g. known_data = shuffle_data(des_mat)), only prepares and skips the shuffling
prepare_data <- function(design_matrix = SS.stack, pred_type = "original", site_no = 1, response = 1, known_data = NULL) {
  
  if (!is.null(known_data)) { # if data is given beforehand, we shall use it instead of the shuffled one
    df_shuffled <- known_data
  } else {
    df_shuffled <- shuffle_data(design_matrix = design_matrix, site_no = site_no, response = response) # function determined above
  }
  
  if (pred_type %in% c("quantile", "normal")) { # turning the predictors into quantiles of their prior distributions

    n_preds <- ncol(df_shuffled) - 1

    idx <- if (n_preds == 5) {
      prior.ind.all
    } else if (n_preds == 10) {
      prior.ind.all_10
    } else if (n_preds == 20) {
      prior.ind.all_20
    } else {
      1:40
    }

    for (i in 1:n_preds) {
      df_shuffled[,i] <- sapply(df_shuffled[,i], function(x) {eval(prior.fn.all$pprior[[idx[i]]], list(q=x))})
    }
  }
  
  if (pred_type == "normal") { # one step further, turning the predictors into standard normal quantile points
    std_norm_values <- apply(df_shuffled[, -ncol(df_shuffled)], MARGIN = c(1,2), function(x) {qnorm(p = x)})
    df_shuffled[, -ncol(df_shuffled)] <- std_norm_values
  }
  
  return(df_shuffled)
  
}


### takes in matrix of experiments, runs new experiments, adds them to that matrix and RETURNS the updated matrix
### (with new experiments)
### CHECK runs.R (can be used to use run_experiment() at background for new experiments)
### FOR DETAILED INFORMATION ABOUT THE FUNCTION & ITS INPUTS!
run_experiments <- function(n_iter, matrix_to_save, knots_list, packages_list, no_params_list, preds_list) {
  
  for (i in 1:n_iter) {
    df_5 <- shuffle_data(des_mat)
    df_10 <- shuffle_data(des_mat_10_10K)
    df_20 <- shuffle_data(des_mat_20_10K)
    df_40 <- shuffle_data(des_mat_40_10K)
    
    for (package in packages_list) {
      for (knots in knots_list) {
        for (pred in preds_list) {
          if (5 %in% no_params_list) {
            metrics_5 <- fit_model(1,1,knots,package,pred,FALSE,0.001,known_data=df_5,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,as.character(knots),metrics_5[[3]],metrics_5[[4]],metrics_5[[1]]+metrics_5[[2]],5,10000)
          }
          if (10 %in% no_params_list) {
            metrics_10 <- fit_model(1,1,knots,package,pred,FALSE,0.001,known_data=df_10,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,as.character(knots),metrics_10[[3]],metrics_10[[4]],metrics_10[[1]]+metrics_10[[2]],10,10000)
          }
          if (20 %in% no_params_list) {
            metrics_20 <- fit_model(1,1,knots,package,pred,FALSE,0.001,known_data=df_20,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,as.character(knots),metrics_20[[3]],metrics_20[[4]],metrics_20[[1]]+metrics_20[[2]],20,10000)
          }
          if (40 %in% no_params_list) {
            metrics_40 <- fit_model(1,1,knots,package,pred,FALSE,0.001,known_data=df_40,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,as.character(knots),metrics_40[[3]],metrics_40[[4]],metrics_40[[1]]+metrics_40[[2]],40,10000)
          }
        }
      }
    }
  }
  
  return(matrix_to_save)
}


### fits a model shape of y = exp(b*x) + c , and RETURNS it
### inputs: y and x as vectors
### used for e.g. to predict rmse by knots or time
exp_decay_offset_mod <- function(y,x) {
  c0 <- 0.5*min(y)
  model0 <- lm(log(y - c0) ~ x)
  a_init <- exp(model0$coefficients[[1]])
  b_init <- model0$coefficients[[2]]
  mod <- nls(y ~ a*exp(b*x)+c, start = list(a = a_init, b = b_init, c = c0), control = nls.control(maxiter = 100))
  return(mod)
}



### for the two following functions, inputs:
### y = observed values
### pred_means = corresponding predictive means
### pred_sds = corresponding predictive standard deviations


### draws a histogram of empirical cdf values of observations y 
### doesn't RETURN anything
### assumes that your model with means pred_means and sd:s pred_sds is right
### based on PIT-theory, these cdf values should be uniformly distributed
### therefore noticeable divergence from Uniform(0,1) suggests that your model isn't right
PIT_histogram <- function(y, pred_means, pred_sds, m = 250) { 
  
  n <- length(y)
  u <- rep(0, n)
  
  pred_sample <- matrix(0, m, n)
  
  for (i in 1:n) {
    pred_sample[,i] <- rnorm(m, pred_means[i], pred_sds[i]) # samples values from predictive marginal univariate normals
    u[i] <- mean(pred_sample[,i] < y[i])
  }
  
  hist(u, breaks = 9, main = "histogram of empirical cdf values at observed values",
       sub = "strong difference from Uniform(0,1) indicates problems in a model")
}


### plots plot.DHARMa figures (QQ-plot, residuals vs predictions) with some p-values too
### doesn't RETURN anything
dharma_figures <- function(y, pred_means, pred_sds, m = 250) { 
  
  n <- length(y)
  fitted <- pred_means
  pred_sample <- matrix(0, n, m)
  
  for (i in 1:n) {
    pred_sample[i,] <- rnorm(m, pred_means[i], pred_sds[i]) # samples values from predictive marginal univariate normals
  }
  
  dharma <- createDHARMa(pred_sample, y, fitted)
  plot(dharma)
}