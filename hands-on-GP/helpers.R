
prepare_and_shuffle_data <- function(design_matrix = SS.stack, pred_type = "original", site_no = 1, response = 1, known_data = NULL) {
  
  if (!is.null(known_data)) { # if data is given beforehand, we shall use it instead of the shuffled one
    df_shuffled <- known_data
  } else {
    df_shuffled <- shuffle_data(design_matrix = design_matrix, site_no = site_no, response = response) # function determined below
  }
  
  if (pred_type %in% c("quantile", "normal")) { # turning the predictors into quantiles of their prior distributions
    df_shuffled$som_respiration_rate <- sapply(df_shuffled$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
    df_shuffled$rdConst <- sapply(df_shuffled$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
    df_shuffled$psnTOpt <- sapply(df_shuffled$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
    df_shuffled$wueConst <- sapply(df_shuffled$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
    df_shuffled$leafGrowth <- sapply(df_shuffled$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})
  }
  
  if (pred_type == "normal") { # one step further, turning the predictors into standard normal quantile points
    std_norm_values <- apply(df_shuffled[, -ncol(df_shuffled)], MARGIN = c(1,2), function(x) {qnorm(p = x)})
    df_shuffled[, -ncol(df_shuffled)] <- std_norm_values
  }
  
  return(df_shuffled)
  
}

shuffle_data <- function(design_matrix = SS.stack, site_no = 1, response = 1) {
  ## check whether your data is just 1 resp 1 site or [[]][[]] type (e.g. SS.stack)?
  if (is.list(design_matrix)) { #SS.stack with 12 sites, 2 responses
    df <- as.data.frame(design_matrix[[site_no]][[response]])
  } else {
    df <- as.data.frame(design_matrix)
  }
  colnames(df)[ncol(df)] <- "likelihood"
  df_shuffled <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ] # shuffles the data
  
  return(df_shuffled)
}
