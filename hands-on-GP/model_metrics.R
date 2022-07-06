library(mlegp)
library(laGP)
library(hetGP)
library(mgcv)
library(tidyverse)
load("whole_stack.Rdata")

##### Fitting GP & predicting ######

model_metrics <- function(site_no = 1, response = 1, n_knots = 100, n_test_points = 150, package = "mlegp",
                          predictor_type = "original", nugget_known = FALSE, nugget = 0.001, gam_type = "bam", data = NULL) {
  
  df <- as.data.frame(SS.stack[[site_no]][[response]])
  colnames(df)[6] <- "likelihood"
  df_shuffled <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  
  if (!is.null(data)) { # if data is given beforehand, we shall use it instead of shuffled
    df_shuffled <- data
  }
  
### prepare quantiles ###  
  if (predictor_type %in% c("quantile", "normal")) {
    df_shuffled$som_respiration_rate <- sapply(df_shuffled$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
    df_shuffled$rdConst <- sapply(df_shuffled$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
    df_shuffled$psnTOpt <- sapply(df_shuffled$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
    df_shuffled$wueConst <- sapply(df_shuffled$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
    df_shuffled$leafGrowth <- sapply(df_shuffled$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})

  }
  
  if (predictor_type == "normal") {
    std_norm_values <- apply(df_shuffled[, -ncol(df_shuffled)], MARGIN = c(1,2), function(x) {qnorm(p = x)})
    df_shuffled[, -ncol(df_shuffled)] <- std_norm_values
  }
  
### set the training & testing sets ##  
  #test_set <- df_shuffled[(n_knots-50):((n_knots-50)+n_test_points), ] #to test how nugget effects when training points also in test set
  #training_set <- df_shuffled[1:n_knots, ]
  test_set <- df_shuffled[(1:n_test_points), ]
  training_set <- df_shuffled[(n_test_points + 1):(n_test_points + n_knots), ]
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
    if (nugget_known) {
      model <- mlegp(X = X, Z = Z, nugget = nugget, nugget.known = 1, verbose = 0)
    } else {
        model <- mlegp(X = X, Z = Z, nugget.known = 0, verbose = 0)
    }
    #print(GPmodel$beta)
  } else if (package == "laGP") {
    d <- darg(list(mle=TRUE, min = sqrt(.Machine$double.eps)), X)
    g <- garg(list(mle=ifelse(nugget_known, FALSE, TRUE)), Z) #mle=FALSE = no estimation for nugget term, just set to 0
    model <- newGPsep(X = X, Z = Z, d = rep(ifelse(predictor_type %in% c("quantile","normal"), 10, d$start), 5), g = ifelse(nugget_known, nugget, 0.01), dK = TRUE) #d$start working badly with quantiles
    if (nugget_known) {
    mle <- mleGPsep(model, param="d", tmin=d$min, tmax=ifelse(predictor_type %in% c("quantile", "normal"), 1e06, d$max), verb = 0)  ###d$max*constant, multiplier to make the upper limit larger so that mle:s dont get stuck to limit? 
    } else {
      #mle <- mleGPsep(model, param="both", tmin=c(d$min, g$min), tmax = c(d$max, 1000), verb = 0)
      mle <- mleGPsep(model, param = "g", tmin=1e-06, tmax=1, verb = 0)
      mleGPsep(model, param = "d", tmin = d$min, tmax=ifelse(predictor_type %in% c("quantile", "normal"), 1e06, d$max), verb = 0)
      g_used <- mle$g
    }
    ### the following to look at the parameters
    #d <- darg(list(mle = TRUE), X) #without setting limits manually
    #print(d)
    #mleGPsep(GPmodel, param = "d", tmin = d$min, tmax = d$max, verb = 1)
    
  } else if (package == "hetGP") {
    if (nugget_known) {
    model <- mleHomGP(X = X, Z = Z, known = list(g = nugget), covtype = "Gaussian")
    } else {
      model <- mleHomGP(X = X, Z = Z, covtype = "Gaussian")
      g_used <- model$g
    }
    #print(GPmodel$theta)
  } else if (package == "mgcv") {
    k <- rep(10, 5)
    if (gam_type == "gam") {
      model <- gam(likelihood ~ s(som_respiration_rate, bs = "cr", k = k[1]) + s(rdConst, bs = "cr", k = k[2]) + s(psnTOpt, bs = "cr", k = k[3]) + s(wueConst, bs = "cr", k = k[4]) + s(leafGrowth, bs = "cr", k = k[5]),
                  data = training_set, method = "REML")
      kcheck <- k.check(model)
      if (min(kcheck[,4]) < 0.01) { ### if k is too low for some predictor, we will double it
        k <- k + k*(kcheck[,4] < 0.01)
        model <- gam(likelihood ~ s(som_respiration_rate, bs = "cr", k = k[1]) + s(rdConst, bs = "cr", k = k[2]) + s(psnTOpt, bs = "cr", k = k[3]) + s(wueConst, bs = "cr", k = k[4]) + s(leafGrowth, bs = "cr", k = k[5]),
                     data = training_set, method = "REML")
        kcheck1 <- kcheck ### comment this off!!!
        kcheck <- k.check(model)
      }

    } else {
        model <- bam(likelihood ~ s(som_respiration_rate, bs = "cr", k = k[1]) + s(rdConst, bs = "cr", k = k[2]) + s(psnTOpt, bs = "cr", k = k[3]) + s(wueConst, bs = "cr", k = k[4]) + s(leafGrowth, bs = "cr", k = k[5]),
                    data = training_set, method = "fREML")
        kcheck <- k.check(model)
        if (min(kcheck[,4]) < 0.01) {
          k <- k + k*(kcheck[,4] < 0.01)
          model <- bam(likelihood ~ s(som_respiration_rate, bs = "cr", k = k[1]) + s(rdConst, bs = "cr", k = k[2]) + s(psnTOpt, bs = "cr", k = k[3]) + s(wueConst, bs = "cr", k = k[4]) + s(leafGrowth, bs = "cr", k = k[5]),
                       data = training_set, method = "fREML")
          kcheck1 <- kcheck
          kcheck <- k.check(model)
        }
      
    }
    # model <- bam() #big data alternative, method = "qREML", bs = "cr" or "ps"
    # print(kcheck)
    # print(gam.check(model))
    # plot.gam(model, pages = 1)
    # if (min(kcheck[,4]) < 0.01) {
    #   print(kcheck1)
    #   print(kcheck)
    # }
  }
  toc <- proc.time()
  time_model <- as.numeric((toc - tic)[3])
  
### predict with the model ###
  tic <- proc.time()
  
  if (package == "mlegp") {
    mod_pred <- predict.gp(model, XX, se.fit = TRUE)
  } else if (package == "laGP") {
    mod_pred <- predGPsep(model, XX, lite = TRUE, nonug = FALSE)
    deleteGPsep(model)
  } else if (package == "hetGP") {
    mod_pred <- predict(model, XX)
  } else if (package == "mgcv") {
    mod_pred <- predict.gam(model, test_set[, -ncol(test_set)], se.fit = TRUE)
  }
  
  toc <- proc.time()
  time_pred <- as.numeric((toc - tic)[3])
  
  if (package %in% c("mlegp", "mgcv")) {
    test_set <- cbind(test_set, mod_pred$fit, mod_pred$se.fit)
    colnames(test_set)[7:ncol(test_set)] <- c("pred", "se")
  } else if (package == "laGP") {
    test_set$pred <- mod_pred$mean
    test_set$se <- sqrt(mod_pred$s2)
    
    #print(mean(sqrt(abs(GP_pred$s2))))
    #print(sqrt(mean((test_set$likelihood - test_set$pred)^2)))
    #print(sum(is.nan(sqrt(GP_pred$s2)))) #why is this producing this much NaN's as GP_pred$s2 aren't NaN?
    #View(cbind(variance = absGP_pred$s2, se = sqrt(GP_pred$s2))) ###########
  } else if (package == "hetGP") {
    test_set$pred <- mod_pred$mean
    test_set$se <- sqrt(mod_pred$sd2 + mod_pred$nugs)
  }
  
### calculate metrics ###  
  if (sum(is.nan(test_set$se)) > 0) {
    nas <- TRUE
  } else {
    nas <- FALSE
  }
  rmse <- sqrt(mean((test_set$likelihood - test_set$pred)^2, na.rm = TRUE)) 
  mean_std_err <- mean((test_set$se), na.rm = TRUE)
  
  metrics <- c("model time" = as.numeric(time_model), "prediction time" = as.numeric(time_pred), "rooted mean square error" = as.numeric(rmse),
               "mean standard error" = as.numeric(mean_std_err), "Nan_SEs" = nas, "sign_p" = ifelse(package == "mgcv", min(kcheck[,4]) < 0.01, 0),
               "est_nug" = ifelse(package %in% c("laGP", "hetGP") & nugget_known == "FALSE", g_used, 0))

  plot(test_set$likelihood, test_set$pred, xlab = "observed", ylab = "predicted", #ylim = c(6000, 13000)
       main = paste0("knots: ", n_knots,", package: ", package, ", pred_type: ", predictor_type))
  abline(0, 1)
  
  
  return(metrics)
  #return(g_used)
  #print(metrics)
  #return(kcheck)
  #print(GPmodel$beta)
  #return(test_set)
  #return(GPmodel)
}

