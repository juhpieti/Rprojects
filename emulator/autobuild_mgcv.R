library(mgcv)

### assume that data is your design matrix and var is the name of your response variabe (e.g. "likelihood" in our case)
### this model will build and return you a gam model with n_interaction interaction terms (default = 0)
### you can also tune your basis (cr = cubic regression as a default)

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
    
    # print(final_formula)
    # print(which.min(unlist(aics[-10])))
    # print(order(unlist(aics), decreasing = FALSE))
    # print(id)
    # print(aics)
    
  } else { # if n_interactions == 0
    final_mod <- if(type=="bam") {
      bam(as.formula(base), data = data, method = "fREML")
    } else {
      gam(as.formula(base), data = data, method = "REML")
    }
  }
  return(final_mod)
}


# 
# tic <- proc.time()
#build_mgcv(training_set, "likelihood", "cr","bam",n_interactions = 3)
# toc <- proc.time()
# print((toc-tic)[3])

#build_mgcv(training_set, "likelihood", "bs", "bam", 4)

# Assuming D_MAT is your design matrix with 6 columns with names of the parameters below and a last column called "llik"
# preds <- c("som_respiration_rate", "psnTOpt", "wueConst", "leafGrowth", "rdConst")
# 
# n  <- length(preds)
# id <- unlist(lapply(1:2, function(i) combn(1:n,i,simplify=FALSE)), recursive=FALSE)
# # also see 
# #id <- unlist(lapply(1:3, function(i) combn(1:n,i,simplify=FALSE)), recursive=FALSE)
# 
# var <- "likelihood"
# # generate formulas
# Formulas <- sapply(id, function(i){
#   if(length(i) == 1){
#     foo <- paste(paste(preds[i],collapse="*"))
#   }else{
#     foo <- paste(paste(preds[i],collapse="+"),"+",paste(paste(preds[i],collapse=":")))
#   }
#   return(paste(var,"~ (", foo, ")"))
# }
# )
# 
# fits <- vector("list", length(Formulas))
# for(i in seq_along(Formulas)){
#   fits[[i]] <- lm(as.formula(Formulas[i]), data=D_MAT)
# }
# #then you can loop over fits and check some properties, is it significant which has the most deviance explained, least penalized etc.
# 
# Formulas

# mod <- build_mgcv(training_set, "likelihood", "cr", "bam", 3)
# gam.check(mod)

