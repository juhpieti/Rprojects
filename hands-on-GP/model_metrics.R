library(mlegp)
library(laGP)
library(hetGP)
library(tidyverse)
#load("toJuho.Rdata")
load("whole_stack.Rdata")

##### Fitting GP & predicting ######

model_metrics <- function(site_no = 1, response = 1, n_knots = 100, n_test_points = 150, package = "mlegp", predictor_type = "normal", data = NULL) {
  df <- as.data.frame(SS.stack[[site_no]][[response]])
  colnames(df)[6] <- "likelihood"
  shuffled_df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  
  if (!is.null(data)) { # if data is given beforehand, we shall use it instead of shuffled
    shuffled_df <- data
  }
  
### prepare quantiles ###  
  if (predictor_type == "quantile") {
    shuffled_df$som_respiration_rate <- sapply(shuffled_df$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
    shuffled_df$rdConst <- sapply(shuffled_df$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
    shuffled_df$psnTOpt <- sapply(shuffled_df$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
    shuffled_df$wueConst <- sapply(shuffled_df$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
    shuffled_df$leafGrowth <- sapply(shuffled_df$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})
  }
  
### set the training & testing sets ##  
  test_set <- shuffled_df[(1:n_test_points), ]
  training_set <- shuffled_df[(n_test_points + 1):(n_test_points + n_knots), ]
  X = training_set[, -ncol(df), drop = FALSE]
  Z = training_set[, ncol(df), drop = ifelse(package %in% c("laGP", "hetGP"), TRUE, FALSE)]
  XX = test_set[ ,-ncol(df), drop = FALSE]
  
  if (package == "hetGP") {
    X <- apply(X, 1:2, function(x) {round(x, digits = 100)}) ### digits fixes the problem with mleHomGP
    XX <- apply(XX, 1:2, function(x) {round(x, digits = 100)}) 
  }
  
### build the model ###  
  tic <- proc.time()
  if (package == "mlegp") {
    GPmodel <- mlegp(X = X, Z = Z, nugget = 0, nugget.known = 1, verbose = 0)
    #print(GPmodel$beta)
  } else if (package == "laGP") {
    d <- darg(list(mle=TRUE, min = sqrt(.Machine$double.eps)), X)
    g <- garg(list(mle=FALSE, start = sqrt(.Machine$double.eps)), Z) #mle=FALSE = no estimation for nugget term, just set to 0
    GPmodel <- newGPsep(X = X, Z = Z, d = rep(d$start, 5), g = g$start, dK = TRUE) #g = 0 or .Machine$double ? g = 0 causes Cholesky decomp error?
    GPmle <- mleGPsep(GPmodel, param="d", tmin=d$min, tmax=d$max)  ###d$max*constant, multiplier to make the upper limit larger so that mle:s dont get stuck to limit? 
    
    ### the following to look at the parameters
    #d <- darg(list(mle = TRUE), X) #without setting limits manually
    #print(d)
    #mleGPsep(GPmodel, param = "d", tmin = d$min, tmax = d$max, verb = 1)
    
  } else if (package == "hetGP") {
    GPmodel <- mleHomGP(X = X, Z = Z, known = list(g = sqrt(.Machine$double.eps)))
    #print(GPmodel$theta)
  }
  toc <- proc.time()
  time_model <- as.numeric((toc - tic)[3])
  
### predict with the model ###
  tic <- proc.time()
  
  if (package == "mlegp") {
    GP_pred <- predict.gp(GPmodel, XX, se.fit = TRUE)
  } else if (package == "laGP") {
    GP_pred <- predGPsep(GPmodel, XX, lite = TRUE, nonug = FALSE)
    deleteGPsep(GPmodel)
  } else if (package == "hetGP") {
    GP_pred <- predict(GPmodel, XX)
  }
  
  toc <- proc.time()
  time_pred <- as.numeric((toc - tic)[3])
  
  if (package == "mlegp") {
    test_set <- cbind(test_set, GP_pred$fit, GP_pred$se.fit)
    colnames(test_set)[7:ncol(test_set)] <- c("pred", "se")
  } else if (package == "laGP") {
    test_set$pred <- GP_pred$mean
    test_set$se <- sqrt(GP_pred$s2)
    
    #print(mean(sqrt(abs(GP_pred$s2))))
    #print(sqrt(mean((test_set$likelihood - test_set$pred)^2)))
    #print(sum(is.nan(sqrt(GP_pred$s2)))) #why is this producing this much NaN's as GP_pred$s2 aren't NaN?
    #View(cbind(variance = absGP_pred$s2, se = sqrt(GP_pred$s2))) ###########
  } else if (package == "hetGP") {
    test_set$pred <- GP_pred$mean
    test_set$se <- sqrt(GP_pred$sd2 + GP_pred$nugs)
  }
  
### calculate metrics ###  
  rmse <- sqrt(mean((test_set$likelihood - test_set$pred)^2, na.rm = TRUE)) 
  mean_std_err <- mean((test_set$se), na.rm = TRUE)
  
  metrics <- c("model time" = as.numeric(time_model), "prediction time" = as.numeric(time_pred), "rooted mean square error" = as.numeric(rmse),
               "mean standard error" = as.numeric(mean_std_err))

  plot(test_set$likelihood, test_set$pred, xlab = "observed", ylab = "predicted", #ylim = c(6000, 13000)
       main = paste0("knots: ", n_knots,", package: ", package, ", pred_type: ", predictor_type))
  abline(0, 1)
  

  #print(metrics)
  #print(GPmodel$beta)
  return(metrics)
  #return(GPmodel)
}

