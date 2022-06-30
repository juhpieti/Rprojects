library(tidyverse)
library(mgcv)
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
load("metrics_matrix_new.Rdata")
load("priorfcn2Juho.Rdata")
source("model_metrics.R")


df <- as.data.frame(SS.stack[[1]][[1]])
colnames(df)[6] <- "likelihood"
shuffled_df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
shuffled_df_01 <- shuffled_df

shuffled_df_01$som_respiration_rate <- sapply(shuffled_df_01$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
shuffled_df_01$rdConst <- sapply(shuffled_df_01$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
shuffled_df_01$psnTOpt <- sapply(shuffled_df_01$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
shuffled_df_01$wueConst <- sapply(shuffled_df_01$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
shuffled_df_01$leafGrowth <- sapply(shuffled_df_01$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})

training_set <- shuffled_df[(151:750), ]
training_set_01 <- shuffled_df_01[(151:250), ]

test_set <- shuffled_df[1:150, ]
test_set_01 <- shuffled_df_01[1:150, ]

X <- training_set[, -ncol(df), drop = FALSE]
X_01 <- training_set_01[, -ncol(df), drop = FALSE]

#X <- apply(X, 1:2, function(x) {round(x, digits = 100)}) ### digits fixes the problem
#X_01 <- apply(X_01, 1:2, function(x) {round(x, digits = 100)})

XX <- test_set[, -ncol(df), drop = FALSE]
#XX <- apply(XX, 1:2, function(x) {round(x, digits = 100)}) 
XX_01 <- test_set_01[, -ncol(df), drop = FALSE]

Z = training_set[, ncol(df), drop = FALSE]
Z_01 = training_set_01[, ncol(df), drop = FALSE]

#reps <- find_reps(X, Z)
tic <- proc.time()
gam_mod <- gam(likelihood ~ s(som_respiration_rate) + s(rdConst) + s(psnTOpt) + s(wueConst) + s(leafGrowth), data = training_set, knots = )
toc <- proc.time()
(toc-tic)[3]
gam_mod
summary(gam_mod)
gam.check(gam_mod)

gam_pred <- predict.gam(gam_mod, test_set[, -ncol(test_set)], se.fit = TRUE)
sqrt(mean((gam_pred$fit - test_set$likelihood)^2))
mean(gam_pred$se.fit)



homGPmodel <- mleHomGP(X = X, Z = Z, known = list(g = sqrt(.Machine$double.eps)), covtype = "Gaussian")#,
#lower = rep(.Machine$double.eps, 5))
homGPpred <- predict(homGPmodel, XX)
toc <- proc.time()
print((tic - toc)[3])

tic <- proc.time()
hetGPmodel <- mleHetGP(X = X, Z = Z, known = list(g = sqrt(.Machine$double.eps)), covtype = "Gaussian")#, settings = list(return.hom = TRUE))
hetGPpred <- predict(hetGPmodel, XX)
toc <- proc.time()
print((tic - toc)[3])


homGPpred <- predict(homGPmodel, XX)
hetGPpred <- predict(hetGPmodel, XX)

homGPpred <- predict(hetGPmodel$modHom, XX)


test_set_mod <- test_set
test_set_mod$pred <- hetGPpred$mean
test_set_mod$se <- sqrt(hetGPpred$sd2 + hetGPpred$nugs)

sqrt(mean((test_set_mod$likelihood - test_set_mod$pred)^2)) #1135
mean(test_set_mod$se) #1194

test_set_mod_01 <- test_set_01
test_set_mod_01$pred = laGPpred_01$mean
test_set_mod_01$se = sqrt(laGPpred_01$s2)

sqrt(mean((test_set_mod_01$likelihood - test_set_mod_01$pred)^2)) #1135
mean(test_set_mod_01$se) #1194

deleteGPsep(laGPmodel)
deleteGPsep(laGPmodel_01)

plot(test_set_mod$likelihood, test_set_mod$pred, ylim = c(6000, 13000))
abline(0, 1)