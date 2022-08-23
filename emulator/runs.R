source("model_metrics.R")

### Run this script in terminal with: "Rscript ~/path/runs.R"
### to produce experiments and save them into a predefined matrix
### it will also print out the time spent to run the experiments

### of course the function run_experiments() can be called in Rstudio to produce experiments into a matrix too (see intro.R)

### HOW TO MODIFY THIS SCRIPT BEFORE RUNNING:
  # 1) load in the matrix / Rdata you want to save your experiments in
  # 2) in the end when you call the function, modify the parameters based on what kind of experiments you want to have
  #    (see function inputs described below)
  # 3) remember to save the matrix with added experiments into a Rdata-file in the end too

### INPUTS:
### n_iter: number of iterations (each iteration has the same data used for different knots, packages etc)
### matrix_to_save: empty matrix or matrix that has the earlier experiments to be expanded with columns:
###                 site_no, response, package, pred_type, no_knots, RMSqE, MStdE, time, no_params, obs   
                 
### knots_list: list of knots you want to experiment with, e.g. knots_list = c(1000, 2000, 3000)
### packages_list: list of packages you want to experiment with, e.g. packages_list = c("laGP, hetGP)
### no_params_list: list of number of parameters you want to experiment with, e.g. no_params_list = c(5, 10, 20)
### preds_list: list of parameter spaces you want to experiment with, e.g. preds_list = c("original", "quantile")


### NOTICEABLE! ###
### this function takes in matrix, runs experiments and adds them to that matrix given as input
### it will then RETURN that modified matrix, so basicly:
### updated_matrix <- run_experiments(n_iter, matrix_before, ...) will just update matrix_before with new experiments and
### save it to updated_matrix


### load in the matrix you want to save your experimentations
load("experiments.Rdata") 

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

  toc <- proc.time()
  print((toc-tic)[[3]])
  
  return(matrix_to_save)
}

### experimenting and saving the experiments in matrix you loaded at the beginning
experiments <- run_experiments(1,experiments,c(5000,6000,7000,8000),c("hetGP"),c(5,10),c("original"))

### saving the updated matrix
save(experiments, file = "experiments.Rdata")

