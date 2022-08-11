
shuffle_data <- function(design_matrix = SS.stack, site_no = 1, response = 1) {
  ## check whether your data is just 1 resp 1 site or [[]][[]] type (e.g. SS.stack)?
  if (length(design_matrix[[1]][[1]]) > 1) { #SS.stack with 12 sites, 2 responses
    df <- as.data.frame(design_matrix[[site_no]][[response]])
  } else {
    df <- as.data.frame(design_matrix)
  }
  colnames(df)[ncol(df)] <- "likelihood"
  df_shuffled <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ] # shuffles the data
  
  return(df_shuffled)
}

prepare_data <- function(design_matrix = SS.stack, pred_type = "original", site_no = 1, response = 1, known_data = NULL) {
  
  if (!is.null(known_data)) { # if data is given beforehand, we shall use it instead of the shuffled one
    df_shuffled <- known_data
  } else {
    df_shuffled <- shuffle_data(design_matrix = design_matrix, site_no = site_no, response = response) # function determined below
  }
  
  if (pred_type %in% c("quantile", "normal")) { # turning the predictors into quantiles of their prior distributions

    n_preds <- ncol(df_shuffled) - 1

    idx <- if (n_preds == 5) {
      prior.ind.all
    } else if (n_preds == 10) {
      prior.ind.all_10
    } else if (n_preds == 20) {
      prior.ind.all_20
    } else {
      1:40
    }

    for (i in 1:n_preds) {
      df_shuffled[,i] <- sapply(df_shuffled[,i], function(x) {eval(prior.fn.all$pprior[[idx[i]]], list(q=x))})
    }
  }
  
  if (pred_type == "normal") { # one step further, turning the predictors into standard normal quantile points
    std_norm_values <- apply(df_shuffled[, -ncol(df_shuffled)], MARGIN = c(1,2), function(x) {qnorm(p = x)})
    df_shuffled[, -ncol(df_shuffled)] <- std_norm_values
  }
  
  return(df_shuffled)
  
}

exp_decay_offset_mod <- function(y,x) { ### fits a model shape of y = exp(b*x) + c
  c0 <- 0.5*min(y)
  model0 <- lm(log(y - c0) ~ x)
  a_init <- exp(model0$coefficients[[1]])
  b_init <- model0$coefficients[[2]]
  mod <- nls(y ~ a*exp(b*x)+c, start = list(a = a_init, b = b_init, c = c0), control = nls.control(maxiter = 100))
  return(mod)
}
