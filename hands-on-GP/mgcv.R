library(tidyverse)
library(mgcv)
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
load("metrics_matrix_new.Rdata")
load("priorfcn2Juho.Rdata")
source("model_metrics.R")


model_metrics(3,1,200,150,"mgcv","quantile",TRUE,1,"gam")

### for loop testatakseen, kuinka paljon p-arvot alittavat esim 0.01? ###


model_metrics()

df <- as.data.frame(SS.stack[[1]][[1]])
colnames(df)[6] <- "likelihood"
shuffled_df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
shuffled_df_01 <- shuffled_df

shuffled_df_01$som_respiration_rate <- sapply(shuffled_df_01$som_respiration_rate, function(x) {eval(prior.fn.all$pprior[[7]], list(q=x))})
shuffled_df_01$rdConst <- sapply(shuffled_df_01$rdConst, function(x) {eval(prior.fn.all$pprior[[10]], list(q=x))})
shuffled_df_01$psnTOpt <- sapply(shuffled_df_01$psnTOpt, function(x) {eval(prior.fn.all$pprior[[25]], list(q=x))})
shuffled_df_01$wueConst <- sapply(shuffled_df_01$wueConst, function(x) {eval(prior.fn.all$pprior[[36]], list(q=x))})
shuffled_df_01$leafGrowth <- sapply(shuffled_df_01$leafGrowth, function(x) {eval(prior.fn.all$pprior[[37]], list(q=x))})

training_set <- shuffled_df[(151:250), ]
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
gam_mod <- gam(likelihood ~ s(som_respiration_rate, bs = "tp", k = 10) + s(rdConst, bs = "tp", k = 10) + s(psnTOpt, bs = "tp", k = 10) + s(wueConst, bs = "tp", k = 20) + s(leafGrowth, bs = "tp", k = 10),
               data = training_set_01, method = "REML")

gam_mod <- gam(likelihood~s(som_respiration_rate, rdConst, psnTOpt, wueConst, leafGrowth), data = training_set_01, method = "REML")



#gam_mod <- bam(likelihood ~ s(som_respiration_rate, bs = "tp") + s(rdConst, bs = "tp") + s(psnTOpt, bs = "tp") + s(wueConst, bs = "tp") + s(leafGrowth, bs = "tp"),
#                data = training_set, method = "fREML", discrete = TRUE)
toc <- proc.time()
(toc-tic)[3]
gam_mod
summary(gam_mod)
gam.check(gam_mod)
k_det <- k.check(gam_mod)

tic <- proc.time()
gam_pred <- predict.gam(gam_mod, test_set_01[, -ncol(test_set_01)], se.fit = TRUE)
toc <- proc.time()
(toc-tic)[3]
sqrt(mean((gam_pred$fit - test_set$likelihood)^2))
mean(gam_pred$se.fit)

plot(gam_pred$fit, test_set$likelihood)
abline(a = 0, b = 1)




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


rmse_normal <- rep(0, 100)
rmse_quantile <- rep(0,100)
mstde_normal <- rep(0,100)
mstde_quantile <- rep(0,100)
time_normal <- rep(0,100)
time_quantile <- rep(0,100)
small_k_normal <- rep(0,100)
small_k_quantile <- rep(0,100)

for (i in 1:100) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_normal <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
                                  package = "mgcv", predictor_type = "normal", data = df)
  metrics_quantile <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
                                    package = "mgcv", predictor_type = "quantile", data = df)
  rmse_normal[i] <- metrics_normal[[3]]
  rmse_quantile[i] <- metrics_quantile[[3]]
  mstde_normal[i] <- metrics_normal[[4]]
  mstde_quantile[i] <- metrics_quantile[[4]]
  time_normal[i] <- metrics_normal[[1]] + metrics_normal[[2]]
  time_quantile[i] <- metrics_quantile[[1]] + metrics_quantile[[2]]
  small_k_normal[i] <- metrics_normal[[6]]
  small_k_quantile[i] <- metrics_quantile[[6]]
}

df_metrics <- data.frame(rmse_normal = rmse_normal, rmse_quantile = rmse_quantile, mstde_normal = mstde_normal,
                         mstde_quantile = mstde_quantile, time_normal = time_normal, time_quantile = time_quantile,
                         "smallK%_normal" = small_k_normal, 'smallK%_quantile' = small_k_quantile)

model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
              package = "mgcv", predictor_type = "normal",gam_type = "bam", data = shuffled_df)

rmse_gam <- rep(0, 100)
rmse_bam <- rep(0,100)
mstde_gam <- rep(0,100)
mstde_bam <- rep(0,100)
time_gam <- rep(0,100)
time_bam <- rep(0,100)
small_k_gam <- rep(0,100)
small_k_bam <- rep(0,100)

for (i in 1:100) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_gam <- model_metrics(site_no = 1, response = 1, n_knots = 600, n_test_points = 150,
                                  package = "mgcv", predictor_type = "quantile", gam_type = "gam", data = df)
  metrics_bam <- model_metrics(site_no = 1, response = 1, n_knots = 600, n_test_points = 150,
                                    package = "mgcv", predictor_type = "quantile", gam_type = "bam", data = df)
  rmse_gam[i] <- metrics_gam[[3]]
  rmse_bam[i] <- metrics_bam[[3]]
  mstde_gam[i] <- metrics_gam[[4]]
  mstde_bam[i] <- metrics_bam[[4]]
  time_gam[i] <- metrics_gam[[1]] + metrics_gam[[2]]
  time_bam[i] <- metrics_bam[[1]] + metrics_bam[[2]]
  small_k_gam[i] <- metrics_gam[[6]]
  small_k_bam[i] <- metrics_bam[[6]]
}

df_metrics <- data.frame(rmse_gam = rmse_gam, rmse_bam = rmse_bam, mstde_gam = mstde_gam,
                         mstde_bam = mstde_bam, time_gam = time_gam, time_bam = time_bam,
                         "smallK%_gam" = small_k_gam, 'smallK%_bam' = small_k_bam)

