library(mgcv)

build_mgcv <- function(data, var, bs, type = "bam", n_interactions = 0) {
  preds <- colnames(data)
  preds <- preds[preds != var]
  n <- length(preds)
  base_k <- rep(10, n)
  base <- paste(paste(var," ~ ", sep = ''), paste('s(', preds,',k=', base_k,",bs='",bs,"'", ')', sep = '', collapse = ' + '), sep = '')
  
  model <- if (type == "bam") {
    bam(as.formula(base), data = data, method = "fREML")
  } else {
    gam(as.formula(base), data = data, method = "REML")
  }
  kcheck <- k.check(model)
  if (min(kcheck[,4]) < 0.05) {
    base_k <- base_k + base_k*(kcheck[,4] < 0.05)
    base <- paste(paste(var," ~ ", sep = ''), paste('s(', preds,',k=', base_k,",bs='",bs,"'", ')', sep = '', collapse = ' + '), sep = '')
  }
  
  ##### interactions? #####
  
  if (n_interactions > 0) {
    id <- combn(1:n,2,simplify=FALSE)
    Formulas <- sapply(id, function(i){
      foo <- paste(preds[i], collapse=", ")
      return(paste(base, paste("ti(", foo, ")", sep = ''), sep = ' + '))
    })
    #print(Formulas)
    
    fits <- vector("list", length = length(Formulas))
    aics <- vector("list", length = length(Formulas))
    for (i in seq_along(Formulas)) {
      fits[[i]] <- if (type == "bam") {
        bam(as.formula(Formulas[i]), data = data, method = "fREML")
      } else {
        gam(as.formula(Formulas[i]), data = data, method = "REML")
      }
      aics[[i]] <- fits[[i]]$aic
    }
    
    interaction_k <- 5
    lowest_aic_idx <- order(unlist(aics), decreasing = FALSE)[1:n_interactions] # just for test I want to include n "best" interaction terms

    interactions <- sapply(id[lowest_aic_idx], function(i){
      foo <- paste(preds[i], collapse=", ")
      return(paste("ti(", foo, ",k=",interaction_k,",bs ='",bs,"'", ")", sep = ""))
    })

    interactions <- paste(interactions, collapse = " + ")
    
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
    
  } else{
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
# hmm <- build_mgcv(training_set, "likelihood", "cr","bam",n_interactions = 0)
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

