library(mgcv)

### assume that data is your design matrix and var is the name of your response variabe (e.g. "likelihood" in our case)
### this function will build and RETURN you a gam model with n_interaction interaction terms (default = 0)
### you can also tune your basis (cr = cubic regression as a default) or set type = "gam" to fit "normal gam"

build_mgcv <- function(data, var, bs = "cr", type = "bam", n_interactions = 0) {
  preds <- colnames(data)
  preds <- preds[preds != var] # extracting the names of just predictors
  n <- length(preds)
  base_k <- rep(10, n) # start with basis dimension k = 10 for every predictor "main effect smooths"
  # generate formula for base model
  base <- paste(paste(var," ~ ", sep = ''), paste('s(', preds,',k=', base_k,",bs='",bs,"'", ')', sep = '', collapse = ' + '), sep = '')
  
  model <- if (type == "bam") {
    bam(as.formula(base), data = data, method = "fREML")
  } else {
    gam(as.formula(base), data = data, method = "REML")
  }
  kcheck <- k.check(model)
  if (min(kcheck[,4]) < 0.05) { # low p-value with edf close to k' (latter one not tested here) may indicate too small of a k'
    base_k <- base_k + base_k*(kcheck[,4] < 0.05) # so we are doubling k in those cases
    # generate new formula for base model with updated basis dimensions
    base <- paste(paste(var," ~ ", sep = ''), paste('s(', preds,',k=', base_k,",bs='",bs,"'", ')', sep = '', collapse = ' + '), sep = '')
  }
  
  ##### interactions #####
  
  if (n_interactions > 0) {
    id <- combn(1:n,2,simplify=FALSE)
    Formulas <- sapply(id, function(i){ # makes Formulas for all base + 1 interaction models
      foo <- paste(preds[i], collapse=", ")
     # print(paste("ti(", foo, ")", sep = '')) ### check this out!
      return(paste(base, paste("ti(", foo, ")", sep = ''), sep = ' + '))
    })
    
    fits <- vector("list", length = length(Formulas))
    aics <- vector("list", length = length(Formulas)) #Akaike information criteria to select the "best" interaction terms
    for (i in seq_along(Formulas)) {
      fits[[i]] <- if (type == "bam") {
        bam(as.formula(Formulas[i]), data = data, method = "fREML")
      } else {
        gam(as.formula(Formulas[i]), data = data, method = "REML")
      }
      aics[[i]] <- fits[[i]]$aic
    }
    
    interaction_k <- 5 # using marginal basis dimensions of 5 for each interaction
    lowest_aic_idx <- order(unlist(aics), decreasing = FALSE)[1:n_interactions] # I want to include n "best" interaction terms

    interactions <- sapply(id[lowest_aic_idx], function(i){ #producting single interaction terms between the n "best"
      foo <- paste(preds[i], collapse=", ")
      return(paste("ti(", foo, ",k=",interaction_k,",bs ='",bs,"'", ")", sep = ""))
    })
    interactions <- paste(interactions, collapse = " + ") # adding them into one formula to add on top of the base model
    
    final_formula <- paste(base, interactions, sep = " + ")
    
    final_mod <- if(type == "bam") {
      bam(as.formula(final_formula), data = data, method = "fREML")
    } else {
      gam(as.formula(final_formula), data = data, method = "REML")
    }
    
    
  } else { # if n_interactions == 0
    final_mod <- if(type=="bam") {
      bam(as.formula(base), data = data, method = "fREML")
    } else {
      gam(as.formula(base), data = data, method = "REML")
    }
  }
  return(final_mod)
}
