library(mlegp)
library(laGP)
library(hetGP)
library(mgcv)
library(tidyverse)

source("helpers.R")

load("whole_stack.Rdata")
load("design_matrix_10K.Rdata")

##### Fitting GP & predicting ######

model_metrics <- function(site_no = 1, response = 1, n_knots = 100, n_test_points = 150, package = "mlegp",
                          predictor_type = "original", nugget_known = FALSE, nugget = 0.001, gam_type = "bam", bs = "tp", data = NULL, design_matrix = SS.stack, verb = 1) {

df <- prepare_and_shuffle_data(design_matrix, predictor_type, site_no, response, data) 

#   df <- as.data.frame(SS.stack[[site_no]][[response]])
#   #df <- as.data.frame(des_mat)
#   colnames(df)[6] <- "likelihood"
#   df_shuffled <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
#   
#   if (!is.null(data)) { # if data is given beforehand, we shall use it instead of the shuffled one
#     df_shuffled <- data
#   }
#   
# ##### prepare predictors #####
#   
#   if (predictor_type %in% c("quantile", "normal")) {
#     df_shuffled$som_respiration_rate <- sapply(df_shuffled$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
#     df_shuffled$rdConst <- sapply(df_shuffled$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
#     df_shuffled$psnTOpt <- sapply(df_shuffled$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
#     df_shuffled$wueConst <- sapply(df_shuffled$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
#     df_shuffled$leafGrowth <- sapply(df_shuffled$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})
# 
#   }
#   
#   if (predictor_type == "normal") {
#     std_norm_values <- apply(df_shuffled[, -ncol(df_shuffled)], MARGIN = c(1,2), function(x) {qnorm(p = x)})
#     df_shuffled[, -ncol(df_shuffled)] <- std_norm_values
#   }
  
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
        if (predictor_type != "original") { ### param = "both" doesn't work with pred_type == "original", maybe due to large range of values within variables
          mle <- mleGPsep(model, param="both", tmin=c(sqrt(.Machine$double.eps), g$min), tmax = c(d$max, 5), verb = 0)
        } else {
            jmleGPsep(model, drange=c(sqrt(.Machine$double.eps), d$max), grange=c(g$min, 5), verb = 0) # this is doing mleGPsep sequentially with
                                                                                                       # param = "d", param == "g" until convergence of MLEs
            # mle <- mleGPsep(model, param = "g", tmin=g$min, tmax=5, verb = 0)
            # mleGPsep(model, param = "d", tmin=sqrt(.Machine$double.eps), tmax=d$max, verb = 0)
        }
      g_used <- ifelse(predictor_type == "original", mle$g, mle$theta[ncol(training_set)])
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
    k <- rep(10, 5) # set the basis dimensions, this (k = 10) is also default for one variable smooth
    if (gam_type == "gam") {
      model <- gam(likelihood ~ s(som_respiration_rate, bs = bs, k = k[1]) + s(rdConst, bs = bs, k = k[2]) + s(psnTOpt, bs = bs, k = k[3]) + s(wueConst, bs = bs, k = k[4]) + s(leafGrowth, bs = bs, k = k[5]),
                  data = training_set, method = "REML")
      kcheck <- k.check(model)
      if (min(kcheck[,4]) < 0.05) { ### if k seems to be too low for some predictor, we will double it
        k <- k + k*(kcheck[,4] < 0.05)
        model <- gam(likelihood ~ s(som_respiration_rate, bs = bs, k = k[1]) + s(rdConst, bs = bs, k = k[2]) + s(psnTOpt, bs = bs, k = k[3]) + s(wueConst, bs = bs, k = k[4]) + s(leafGrowth, bs = bs, k = k[5]),
                     data = training_set, method = "REML")
        kcheck1 <- kcheck ### comment this off!!!
        kcheck <- k.check(model)
      }

    } else { # gam_type == "bam" (we decided to use this as a default)
             # with bam one should use method = "fREML" and bs = e.g. "cr", "ps" (avoid default bs = "tp")
        model <- bam(likelihood ~ s(som_respiration_rate, bs = bs, k = k[1]) + s(rdConst, bs = bs, k = k[2]) + s(psnTOpt, bs = bs, k = k[3]) + s(wueConst, bs = bs, k = k[4]) + s(leafGrowth, bs = bs, k = k[5]),
                    data = training_set, method = "fREML")
        kcheck <- k.check(model)
        if (min(kcheck[,4]) < 0.05) {
          k <- k + k*(kcheck[,4] < 0.05)
          model <- bam(likelihood ~ s(som_respiration_rate, bs = bs, k = k[1]) + s(rdConst, bs = bs, k = k[2]) + s(psnTOpt, bs = bs, k = k[3]) + s(wueConst, bs = bs, k = k[4]) + s(leafGrowth, bs = bs, k = k[5]),
                       data = training_set, method = "fREML")
          kcheck1 <- kcheck
          kcheck <- k.check(model)
        }
        
        pred_names <- colnames(training_set)[1:5]
        if (min(kcheck[,4]) < 0.05) {
          small_k_names <- pred_names[kcheck[,4] < 0.05]
        } else {
          small_k_names <- "None"
        }
        # # print(pred_names)
        # print(min(kcheck[,4]) < 0.05)
        # print(small_k_names)
        #print(kcheck)
    }
     #print(kcheck)
    #print(gam.check(model))
     #plot.gam(model, pages = 1)
    # if (min(kcheck[,4]) < 0.05) {
    #   print(paste0(predictor_type, ", ",bs))
    #   print(kcheck1)
    #   print(kcheck)
    # }
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
    colnames(test_set)[7:ncol(test_set)] <- c("pred", "se")
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
                 "mean standard error" = as.numeric(mean_std_err), "Nan_SEs" = nas, "sign_p" = ifelse(package == "mgcv", min(kcheck[,4]) < 0.05, FALSE),
                 "est_nug" = ifelse(package %in% c("laGP", "hetGP") & nugget_known == "FALSE", g_used, 0),
                 "failed_predictor" = ifelse(package == "mgcv", small_k_names, 0))
  } else {
    metrics <- c("model time" = as.numeric(time_model), "prediction time" = as.numeric(time_pred), "rooted mean square error" = as.numeric(rmse),
                 "mean standard error" = as.numeric(mean_std_err))
  }
  
  plot(test_set$likelihood, test_set$pred, xlab = "observed", ylab = "predicted",
       main = paste0("site: ",site_no,", resp: ",response,", knots: ", n_knots,", package: ", package, ", pred_type: ", predictor_type))
  abline(0, 1)
  
  
  return(metrics)
  #return(g_used)
  #print(metrics)
  #return(kcheck)
  #return(model)
  #return(mle)
}

