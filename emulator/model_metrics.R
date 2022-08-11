# loading in packages used
library(mlegp)
library(laGP)
library(hetGP)
library(mgcv)

# laoding in functions used inside this fit_model function
source("helpers.R")
source("autobuild_mgcv.R")

# loading in the design matrices
load("data/whole_stack.Rdata")
load("data/design_matrix_10K.Rdata")
load("data/design_matrix_10x10K.Rdata")
load("data/design_matrix_20x10K.Rdata")
load("data/design_matrix_40x10K.Rdata")

# loading the prior distributions & indexes
load("data/priorfcn2Juho.Rdata")
load("data/prior.ind.all_10.Rdata")
load("data/prior.ind.all_20.Rdata")


##### Fitting GP & predicting ######

fit_model <- function(site_no = 1, response = 1, n_knots = 100, n_test_points = 150, package = "mlegp",
                          pred_type = "original", nugget_known = FALSE, nugget = 0.001, gam_type = "bam", bs = "cr",
                          n_interactions = 0, design_matrix = SS.stack, known_data = NULL, verb = 0) {
  
  df <- prepare_and_shuffle_data(design_matrix, pred_type, site_no, response, known_data)  #function from helpers.R
  # NOTE: if known data given (known_data != NULL), function isn't shuffling, can be used for comparisons

##### set the training & testing sets #####
  test_set <- df[(1:n_test_points), ]
  training_set <- df[(n_test_points + 1):(n_test_points + n_knots), ]
  
  X = training_set[, -ncol(df), drop = FALSE] # predictor matrix for training
  Z = training_set[, ncol(df), drop = ifelse(package %in% c("laGP", "hetGP"), TRUE, FALSE)] # vector of observations for training
  XX = test_set[ ,-ncol(df), drop = FALSE] # predictor matrix for testing
  
  if (package == "hetGP") {
    X <- apply(X, 1:2, function(x) {round(x, digits = 100)}) # rounding (with a large number of digits)digits
    XX <- apply(XX, 1:2, function(x) {round(x, digits = 100)}) # fixed a problem "length(Z) should be equal to sum(mult)" with mleHomGP
                                                               # where sum(mult) = nrow(X) in our case with no replicates 
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
      
      ### darg and garg helps you to automate model building by giving you parameters to use
      ### with mleGP functions (e.g. ranges to search the length-scale or nugget parameters from, initial values)
      ### d refers to length-scale parameters, g refers to nugget 
      d <- darg(list(mle=TRUE), X)
      g <- garg(list(mle=ifelse(nugget_known, FALSE, TRUE)), Z) # mle=FALSE = no estimation for nugget term, set to starting value
      
      model <- newGPsep(X = X, Z = Z, d = rep(d$start, ncol(X)), g = ifelse(nugget_known, nugget, 0.1), dK = TRUE) 
      
      if (nugget_known) {
        mle <- mleGPsep(model, param="d", tmin=sqrt(.Machine$double.eps), tmax=d$max, verb = 0) 
      } else {
        if (pred_type != "original") { # param = "both" doesn't work with pred_type == "original", maybe due to large range of values within variables
          mle <- mleGPsep(model, param="both", tmin=c(sqrt(.Machine$double.eps), g$min), tmax = c(d$max, 5), verb = 0)
        } else {
            # mle <- jmleGPsep.R(model, N = 2, drange=c(sqrt(.Machine$double.eps), d$max), grange=c(g$min, 5), verb = 1)
            mle <- jmleGPsep(model, drange=c(sqrt(.Machine$double.eps), d$max), grange=c(g$min, 5), verb = 0) # this is doing mleGPsep sequentially with
                                                                                                              # param = "d", param == "g" until convergence of MLEs
            ### running the following would estimate wrt. both just once and leave it there,
            ### saving time with the cost of losing accuracy
            
            # mleGPsep(model, param = "d", tmin=sqrt(.Machine$double.eps), tmax=d$max, verb = 0) 
            # mle <- mleGPsep(model, param = "g", tmin=g$min, tmax=5, verb = 0)
        }
      g_used <- ifelse(pred_type == "original", mle$g, mle$theta[ncol(training_set)]) # if jmleGPsep used
      #g_used <- ifelse(pred_type == "original", mle$mle[length(mle$mle)-2], mle$theta[ncol(training_set)]) # if jmleGPsep.R
      }
      # print(d) # for checking how did the estimation work out, use with verb = 1 to see the estimated values
                 # sometimes it happens that the estimated values just end up at max$d, or stays at the starting value, which is suspicious
       
  ### build the hetGP-model ###
  } else if (package == "hetGP") {
    if (nugget_known) {
      model <- mleHomGP(X = X, Z = Z, known = list(g = nugget), covtype = "Gaussian")
    } else {
      model <- mleHomGP(X = X, Z = Z, covtype = "Gaussian")
      g_used <- model$g
    }
    #print(model$nit_opt) # how many iterations to get the mle-values for parameters
    
    #print(model$theta) # to compare with other models
    
  ### build the mgcv-model ###
  } else if (package == "mgcv") {
    
    model <- build_mgcv(training_set, "likelihood", bs = bs, type = gam_type, n_interactions = n_interactions) #function from autobuild_mgcv.R
    # model <- build_mgcv(training_set, colnames(X)[ncol(X)], bs = bs, type = gam_type, n_interactions = n_interactions)
  }
  
  toc <- proc.time()
  time_model <- as.numeric((toc - tic)[3])
  
##### predict with the model #####
  
  tic <- proc.time()
  
  if (package == "mlegp") {
    mod_pred <- predict.gp(model, XX, se.fit = TRUE)
  } else if (package == "laGP") {
    mod_pred <- predGPsep(model, XX, lite = TRUE)
    deleteGPsep(model) # deletes GP object from C-side: "Any function calling newGP or newGPsep will require destruction
                       # via these functions or there will be a memory leak"s
  } else if (package == "hetGP") {
    mod_pred <- predict(model, XX)
  } else if (package == "mgcv") {
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
    test_set$se <- sqrt(mod_pred$sd2 + mod_pred$nugs) # "The full predictive variance corresponds to the sum of sd2 and nugs."
                                                      # - Documentation of predict.hetGP
                                                      #  The output from predict separates variance in terms of mean (p$sd2)
                                                      #  and residual (nugget-based p$nugs) estimates, which we combine
                                                      #  to get the full predictive uncertainty of the response Y (x) | DN .
                                                      #  - Vignette of hetGP
  
    # print(scores(model, XX, test_set$likelihood, return.rmse = TRUE))
  }
  
##### calculate the metrics #####
  
  if (sum(is.nan(test_set$se)) > 0) { # there might be NaNs when using nugget too close to zero
                                      # (variances went negative, don't know the exact reason)
                                      # this nas variable was used to notice these problematic cases
    nas <- TRUE
  } else {
    nas <- FALSE
  }
  
  rmse <- sqrt(mean((test_set$likelihood - test_set$pred)^2, na.rm = TRUE)) 
  mstde <- mean((test_set$se), na.rm = TRUE)
  
  if (verb == 1) {
    metrics <- list("mod_time" = as.numeric(time_model), "pred_time" = as.numeric(time_pred), "rmse" = as.numeric(rmse),
                 "mstde" = as.numeric(mstde), "observed" = test_set$likelihood, "pred" = test_set$pred, "pred_se" = test_set$se, "Nan_SEs" = nas,
                 "est_nug" = ifelse(package %in% c("laGP", "hetGP") & nugget_known == "FALSE", g_used, 0))
  } else {
    metrics <- c("mod_time" = as.numeric(time_model), "pred_time" = as.numeric(time_pred), "rmse" = as.numeric(rmse),
                 "mstde" = as.numeric(mstde))
  }
  
  # plot(test_set$likelihood, test_set$pred, xlab = "observed", ylab = "predicted",
  #      main = paste0("site: ",site_no,", resp: ",response,", knots: ", n_knots,", package: ", package, ", pred_type: ", pred_type))
  plot(test_set$pred, test_set$likelihood, xlab = "predicted", ylab = "observed",
       main = paste0("package: ",package,", knots: ",n_knots,", predictors: ",ncol(X),", type: ", pred_type),
       sub = paste0("data: ",nrow(design_matrix)," observations of response ", response, " from site ",site_no))
  abline(0, 1)
  
  return(metrics)
  #return(model)
}

