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

### INPUTS of run_experiment():
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


### load in the matrix you want to save your experiments
load("data/experiments.Rdata") 

### experimenting and saving the experiments in matrix you loaded at the beginning
tic <- proc.time() # to see how much time this iteration set took

experiments <- run_experiments(1,experiments,c(5000,6000,7000,8000),c("hetGP"),c(5,10),c("original")) # function from helpers.R
toc <- proc.time()

print((toc-tic)[[3]])

### saving the updated matrix
save(experiments, file = "data/experiments.Rdata")

