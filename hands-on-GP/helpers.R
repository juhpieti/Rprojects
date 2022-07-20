prepare_and_shuffle_data <- function(des_matrix, pred_type, site_no, response, data = NULL) {
  ## check whether your data is just 1 resp 1 site or [[]][[]] type?
  if (is.list(des_matrix)) { #SS.stack with 12 sites, 2 responses
    print("hep")
    df <- as.data.frame(des_matrix[[site_no]][[response]])
  } else {
    df <- as.data.frame(des_matrix)
  }
  colnames(df)[6] <- "likelihood"
  df_shuffled <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  
  if (!is.null(data)) { # if data is given beforehand, we shall use it instead of the shuffled one
    df_shuffled <- data
  }
  
  if (pred_type %in% c("quantile", "normal")) {
    df_shuffled$som_respiration_rate <- sapply(df_shuffled$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
    df_shuffled$rdConst <- sapply(df_shuffled$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
    df_shuffled$psnTOpt <- sapply(df_shuffled$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
    df_shuffled$wueConst <- sapply(df_shuffled$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
    df_shuffled$leafGrowth <- sapply(df_shuffled$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})
    
  }
  
  if (pred_type == "normal") {
    std_norm_values <- apply(df_shuffled[, -ncol(df_shuffled)], MARGIN = c(1,2), function(x) {qnorm(p = x)})
    df_shuffled[, -ncol(df_shuffled)] <- std_norm_values
  }
  
  return(df_shuffled)
}