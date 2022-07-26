library(mlegp)
library(laGP)
library(hetGP)
library(mgcv)
library(tidyverse)

source("helpers.R")
source("autobuild_mgcv.R")

load("whole_stack.Rdata")
load("design_matrix_10K.Rdata")
load("design_matrix_10x10K.Rdata")

##### Fitting GP & predicting ######

model_metrics <- function(site_no = 1, response = 1, n_knots = 100, n_test_points = 150, package = "mlegp",
                          pred_type = "original", nugget_known = FALSE, nugget = 0.001, gam_type = "bam", bs = "cr",
                          n_interactions = 0, design_matrix = SS.stack, known_data = NULL, verb = 1) {
  
  df <- prepare_and_shuffle_data(design_matrix, pred_type, site_no, response, known_data)  #function from Helpers.R

##### set the training & testing sets #####
  
  test_set <- df[(1:n_test_points), ]
  training_set <- df[(n_test_points + 1):(n_test_points + n_knots), ]
  
  X = training_set[, -ncol(df), drop = FALSE]
  Z = training_set[, ncol(df), drop = ifelse(package %in% c("laGP", "hetGP"), TRUE, FALSE)]
  XX = test_set[ ,-ncol(df), drop = FALSE]
  
  if (package == "hetGP") {
    X <- apply(X, 1:2, function(x) {round(x, digits = 100)}) ### digits fixes the problem with mleHomGP
    XX <- apply(XX, 1:2, function(x) {round(x, digits = 100)}) 
  }
  
##### build the model #####
  
  tic <- proc.time()
  
  ### build the mlegp-model ###
  if (package == "mlegp") {
    if (nugget_known) {
      model <- mlegp(X = X, Z = Z, nugget = nugget, nugget.known = 1, verbose = 0)
    } else {
      model <- mlegp(X = X, Z = Z, nugget.known = 0, verbose = 0)
      g_used <- model$nugget
    }
    #print(model$beta) #to compare between other packages
    
  ### build the laGP-model ###
  } else if (package == "laGP") {
      d <- darg(list(mle=TRUE), X)
      g <- garg(list(mle=ifelse(nugget_known, FALSE, TRUE)), Z) #mle=FALSE = no estimation for nugget term, set to starting value
      model <- newGPsep(X = X, Z = Z, d = rep(d$start, 5), g = ifelse(nugget_known, nugget, 0.1), dK = TRUE) 
      if (nugget_known) {
        mle <- mleGPsep(model, param="d", tmin=sqrt(.Machine$double.eps), tmax=d$max, verb = 0) 
      } else {
        if (pred_type != "original") { ### param = "both" doesn't work with pred_type == "original", maybe due to large range of values within variables
          mle <- mleGPsep(model, param="both", tmin=c(sqrt(.Machine$double.eps), g$min), tmax = c(d$max, 5), verb = 0)
        } else {
            #mle <- jmleGPsep(model, drange=c(sqrt(.Machine$double.eps), d$max), grange=c(g$min, 5), verb = 0) # this is doing mleGPsep sequentially with
                                                                                                       # param = "d", param == "g" until convergence of MLEs
            mle <- mleGPsep(model, param = "g", tmin=g$min, tmax=5, verb = 0)
            mleGPsep(model, param = "d", tmin=sqrt(.Machine$double.eps), tmax=d$max, verb = 0)
        }
      g_used <- ifelse(pred_type == "original", mle$g, mle$theta[ncol(training_set)])
      }
      #print(d)
       
  ### build the hetGP-model ###
  } else if (package == "hetGP") {
    if (nugget_known) {
      model <- mleHomGP(X = X, Z = Z, known = list(g = nugget), covtype = "Gaussian")
    } else {
      model <- mleHomGP(X = X, Z = Z, covtype = "Gaussian")
      g_used <- model$g
    }
    #print(model$nit_opt)
    
    #print(model$theta) # to compare with other models
    
  ### build the mgcv-model ###
  } else if (package == "mgcv") {
    
    model <- build_mgcv(training_set, "likelihood", bs = bs, type = gam_type, n_interactions = n_interactions) #function from autobuild_mgcv.R
    
  }
  
  toc <- proc.time()
  time_model <- as.numeric((toc - tic)[3])
  
##### predict with the model #####
  
  tic <- proc.time()
  
  if (package == "mlegp") {
    mod_pred <- predict.gp(model, XX, se.fit = TRUE)
  } else if (package == "laGP") {
    mod_pred <- predGPsep(model, XX, lite = TRUE, nonug = FALSE)
    deleteGPsep(model)
  } else if (package == "hetGP") {
    mod_pred <- predict(model, XX)
  } else if (package == "mgcv") {
    #mod_pred <- predict.gam(model, test_set[, -ncol(test_set)], se.fit = TRUE)
    mod_pred <- predict.gam(model, XX, se.fit = TRUE)
  }
  
  toc <- proc.time()
  time_pred <- as.numeric((toc - tic)[3])
  
  if (package %in% c("mlegp", "mgcv")) {
    test_set <- cbind(test_set, mod_pred$fit, mod_pred$se.fit)
    colnames(test_set)[(ncol(test_set) - 1):ncol(test_set)] <- c("pred", "se")
  } else if (package == "laGP") {
    test_set$pred <- mod_pred$mean
    test_set$se <- sqrt(mod_pred$s2)
  } else if (package == "hetGP") {
    test_set$pred <- mod_pred$mean
    test_set$se <- sqrt(mod_pred$sd2 + mod_pred$nugs)
  }
  
##### calculate the metrics #####
  
  if (sum(is.nan(test_set$se)) > 0) { ### there might be NaNs when using nugget too close to zero (don't know the exact reason)
    nas <- TRUE
  } else {
    nas <- FALSE
  }
  
  rmse <- sqrt(mean((test_set$likelihood - test_set$pred)^2, na.rm = TRUE)) 
  mean_std_err <- mean((test_set$se), na.rm = TRUE)
  
  if (verb == 1) {
    metrics <- list("model time" = as.numeric(time_model), "prediction time" = as.numeric(time_pred), "rooted mean square error" = as.numeric(rmse),
                 "mean standard error" = as.numeric(mean_std_err), "Nan_SEs" = nas,
                 "est_nug" = ifelse(package %in% c("laGP", "hetGP") & nugget_known == "FALSE", g_used, 0))
  } else {
    metrics <- c("model time" = as.numeric(time_model), "prediction time" = as.numeric(time_pred), "rooted mean square error" = as.numeric(rmse),
                 "mean standard error" = as.numeric(mean_std_err))
  }
  
  plot(test_set$likelihood, test_set$pred, xlab = "observed", ylab = "predicted",
       main = paste0("site: ",site_no,", resp: ",response,", knots: ", n_knots,", package: ", package, ", pred_type: ", pred_type))
  abline(0, 1)
  
  
  return(metrics)
  #return(model)
}

