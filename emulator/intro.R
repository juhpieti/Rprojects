source("model_metrics.R")
source("plotting_functions.R")

### the main function is fit_model() from model_metrics.R
### basicly fit_model() fits you a chosen model to chosen design matrix and returns (as well as prints)
### basic model metrics: time consumptions (fit & predict), rooted mean square error (accuracy), mean standard error (certainty)
### it also plots you observed vs predicted as a default, other diagnostics plots with option "diagnostics = TRUE"
### look up that R-file (model_metrics.R) for more detailed description of fit_model() and its inputs!



##### BASIC USE OF FIT_MODEL() #####

### let's see an example:
fit_model(site_no = 1, response = 1, n_knots = 100, package = "mlegp", pred_type = "original",
          nugget_known = TRUE, nugget = 0, design_matrix = SS.stack, diagnostics = FALSE, verb = 0)

### that call fitted a GP (with mlegp) to design matrix from SS.stack (site 1, response 1) with 100 knots,
### calculated some metrics, printed them and plotted observed values vs predicted values with 45 degree line

### another one:
set.seed(12345)

output <- fit_model(n_knots = 1000, package = "hetGP", pred_type = "quantile", 
          nugget_known = FALSE, design_matrix = des_mat_10_10K, diagnostics = TRUE, verb = 1)

output$est_nug # 0.01437

### that call fitted a GP (with hetGP) using des_mat_10_10K (10 000 x 11 matrix) with 1000 knots
### no need to think about site_no and response as larger data are just from site 1 with response 1
### we also changed from original parameter space to quantile space and used estimated nugget instead of 0
### diagnostics = TRUE gave us a couple extra plots to examine the suitability of the model fitted
### verb = 1 gave us a bit more information about the fit & performance, e.g. estimated nugget = 0.01625



##### SHUFFLE_DATA() TO PRODUCE COMPARABLE RUNS #####

### I'll also demonstrate use of shuffle_data() to produce comparable results between packages
### fit_model uses this function to shuffle the design matrix (then it picks 20% first to test, then n_knots following to train)
### as a default.
### But we can bypass that by shuffling data ourselves first, then giving it as input to fit_model():
df <- shuffle_data(des_mat)
head(df) # rows in random order
fit_model(n_knots = 600, package = "laGP", pred_type = "normal", nugget_known = FALSE, known_data = df)

### now as we run again using the same "known data" (meaning: re-shuffled), we get the same results
fit_model(n_knots = 600, package = "laGP", pred_type = "normal", nugget_known = FALSE, known_data = df)

### then we can check how for example hetGP is doing with the exactly same settings
fit_model(n_knots = 600, package = "hetGP", pred_type = "normal", nugget_known = FALSE, known_data = df)



##### RUN_EXPERIMENTS TO RUN & SAVE EXPERIMENTS #####

### the following is now probably a bit messy / silly way of doing things, but here's what I end up with for experimentation:

### first create a data frame to save experiment results in (these columns are expected from data frame given
### as input for my plotting functions, we will have one example below)

df_experiments <- data.frame(site_no = NA, response = NA, package = NA, pred_type = NA, no_knots = NA,
                             RMSqE = NA, MStdE = NA, time = NA, no_params = NA, obs = NA)
(df_experiments <- df_experiments[-1, ]) # delete the first row for just empty data frame with column names

df_experiments <- run_experiments(n_iter = 3, matrix_to_save = df_experiments, knots_list = c(100,200,300,400,500,600),
                                  packages_list = c("hetGP", "laGP", "mgcv"), no_params_list = c(5),
                                  preds_list = c("quantile", "normal"))

### run_experiments() does n_iter iterations of model runs which settings given as input
### e.g. here we made 12 runs (6 different number of knots, 2 different param spaces) for each package (hetGP, laGP, mgcv)
### so in total 12*3 = 36 runs
### n_iter = 3 told to do this process 3 times, giving 3*36 = 108 runs
### no_params_list = c(5) told to experiment only 5-parameter-data. We could give e.g c(5,10,20,40) to explore dimensions.
### !check runs.R script (used to run experiments on background) for more information about run_experiments()!

nrow(df_experiments) # 108 as expected
head(df_experiments, 10)




##### VISUALIZING THE RESULTS OF EXPERIMENTS #####


### then we can use some functions from plotting_functions.R to visualize the experiments
compare_packages_plot(matrix = df_experiments, pred = "quantile", metric = "RMSqE", n_obs = 10000, n_params = 5)
compare_param_spaces_plot(matrix = df_experiments, n_obs = 10000)

### If we had experiments with different amount of parameters too, we could use compare_no_params_plot()
### to visualize differences between those


##### END #####
