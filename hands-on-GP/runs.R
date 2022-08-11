### run this script in terminal with: "Rscript ~/FMI/Rprojects/hands-on-GP/helloworld.R"

source("model_metrics.R")
load("comp_mat_hetGP_dim_2k.Rdata")
load("laGP_mat.Rdata")


run_experiments <- function(n_iter, matrix_to_save, knots_list, packages_list, no_params_list, preds_list) {
  
  tic <- proc.time() ### to see how much time this iteration set took

  for (i in 1:n_iter) {
    df_5 <- shuffle_data(des_mat)
    df_10 <- shuffle_data(des_mat_10_10K)
    df_20 <- shuffle_data(des_mat_20_10K)
    df_40 <- shuffle_data(des_mat_40_10K)

    for (package in packages_list) {
      for (knots in knots_list) {
        for (pred in preds_list) {
          if (5 %in% no_params_list) {
            metrics_5 <- fit_model(1,1,knots,2000,package=package,pred_type=pred,FALSE,0.001,known_data=df_5,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,knots,metrics_5[[3]],metrics_5[[4]],metrics_5[[1]]+metrics_5[[2]],5,10000)
          }
          if (10 %in% no_params_list) {
            metrics_10 <- fit_model(1,1,knots,2000,package=package,pred_type=pred,FALSE,0.001,known_data=df_10,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,knots,metrics_10[[3]],metrics_10[[4]],metrics_10[[1]]+metrics_10[[2]],10,10000)
          }
          if (20 %in% no_params_list) {
            metrics_20 <- fit_model(1,1,knots,2000,package=package,pred_type=pred,FALSE,0.001,known_data=df_20,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,knots,metrics_20[[3]],metrics_20[[4]],metrics_20[[1]]+metrics_20[[2]],20,10000)
          }
          if (40 %in% no_params_list) {
            metrics_40 <- fit_model(1,1,knots,2000,package=package,pred_type=pred,FALSE,0.001,known_data=df_40,verb=0)
            matrix_to_save[nrow(matrix_to_save)+1, ] <- list(1,1,package,pred,knots,metrics_40[[3]],metrics_40[[4]],metrics_40[[1]]+metrics_40[[2]],40,10000)
          }
        }
      }
    }
  }

  toc <- proc.time()
  print((toc-tic)[[3]])
  
  return(matrix_to_save)
  
}

#comp_mat_hetGP_2k <- run_experiments(1, comp_mat_hetGP_2k)
#save(comp_mat_hetGP_2k, file = "comp_mat_hetGP_dim_2k.Rdata")

laGP_mat <- run_experiments(2,laGP_mat,c(1500,2000),c("laGP"),c(5,10,20),c("original"))
save(laGP_mat, file = "laGP_mat.Rdata")

