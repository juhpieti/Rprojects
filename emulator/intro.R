source("model_metrics.R")
source("runs.R")
source("plotting_functions.R")

### the main function is fit_model() from model_metrics.R
### basicly fit_model() fits you a chosen model to chosen design matrix and returns (as well as prints)
### basic model metrics: time consumption, rooted mean square error (accuracy), mean standard error (certainty)
### it also plots you observed vs predicted as a default, other diagnostics plots with option "diagnostics = TRUE"
### look up that R-file for more detailed description of fit_model() and its inputs!

### let's see an example:
fit_model(site_no = 1, response = 1, n_knots = 100, package = "mlegp", pred_type = "original",
          nugget_known = TRUE, nugget = 0, design_matrix = SS.stack, diagnostics = FALSE, verb = 0)

### that call fitted a GP (with mlegp) to design matrix from SS.stack (site 1, response 1) with 100 knots,
### calculated some metrics, printed them and plotted observed values vs predicted values with 45 degree line

### another one:
fit_model(n_knots = 1000, package = "hetGP", pred_type = "quantile", 
          nugget_known = FALSE, design_matrix = des_mat_10_10K, diagnostics = TRUE, verb = 1)

### that call fitted a GP (with hetGP) to des_mat_10_10K (10 000 x 11 matrix) with 1000 knots
### no need to think about site_no and response as larger data are just from site 1 with response 1
### we also changed from original parameter space to quantile space and used estimated nugget
### diagnostics = TRUE gave us a couple extra plots to examine the suitability of the model fitted
### verb = 1 gave us a bit more information about the fit & performance, e.g. estimated nugget = 0.0099

###### matrix of experiments? use runs.R to create experiments? use plot functions to give plots?