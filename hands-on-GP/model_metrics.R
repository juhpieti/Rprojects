library(mlegp)
library(tidyverse)
load("toJuho.Rdata")
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")

##### Fitting GP, predicting ######

model_metrics <- function(site_no = 1, response = 1, n_knots = 100, n_test_points = 150) {
  df <- as.data.frame(SS.stack[[site_no]][[response]])
  colnames(df)[6] <- "likelihood"
  shuffled_df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  test_set <- shuffled_df[(1:n_test_points), ]
  training_set <- shuffled_df[(n_test_points + 1):(n_test_points + n_knots), ]
  
  tic <- proc.time()
  GPmodel <- mlegp(X = training_set[, -ncol(df), drop = FALSE], Z = training_set[, ncol(df), drop = FALSE], 
                   nugget = 0, nugget.known = 1, verbose = 0, )
  toc <- proc.time()
  time_model = as.numeric((toc - tic)[3])
  
  tic <- proc.time()
  GP_pred = predict.gp(GPmodel, test_set[ ,-ncol(df), drop = FALSE], se.fit = TRUE)
  
  toc <- proc.time()
  time_pred <- as.numeric((toc - tic)[3])
  
  
  test_set <- cbind(test_set, GP_pred$fit, GP_pred$se.fit)
  colnames(test_set)[7:ncol(test_set)] <- c("pred", "se")
  
  rmse <- sqrt(mean((test_set$likelihood - test_set$pred)^2)) 
  mean_std_err <- mean((test_set$se))
  
  metrics <- c("model time" = time_model, "prediction time" = time_pred, "rooted mean square error" = rmse,
               "mean standard error" = mean_std_err)
  
  plot(GPmodel)
  return(metrics)
  
}

