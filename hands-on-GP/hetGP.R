library(tidyverse)
library(hetGP)
load("whole_stack.Rdata")
load("metrics_matrix.Rdata")
load("metrics_matrix_new.Rdata")
load("priorfcn2Juho.Rdata")
source("model_metrics.R")

lol <- model_metrics(1,1,100,150,"mlegp","normal")

model_metrics(1,1,300,150,"hetGP","normal",TRUE,sqrt(.Machine$double.eps),shuffled_df)
model_metrics(1,1,300,150,"hetGP","normal",TRUE,0.1,shuffled_df)
model_metrics(1,1,300,150,"hetGP","normal",TRUE,0.01,shuffled_df)
model_metrics(1,1,300,150,"hetGP","normal",TRUE,0.001,shuffled_df)
model_metrics(1,1,300,150,"hetGP","normal",FALSE,0.1,shuffled_df)





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

X <- apply(X, 1:2, function(x) {round(x, digits = 100)}) ### digits fixes the problem
X_01 <- apply(X_01, 1:2, function(x) {round(x, digits = 100)})

XX <- test_set[, -ncol(df), drop = FALSE]
XX <- apply(XX, 1:2, function(x) {round(x, digits = 100)}) 
XX_01 <- test_set_01[, -ncol(df), drop = FALSE]


Z = training_set[, ncol(df), drop = TRUE]
Z_01 = training_set_01[, ncol(df), drop = TRUE]

#reps <- find_reps(X, Z)
tic <- proc.time()
homGPmodel <- mleHomGP(X = X, Z = Z, known = list(g = sqrt(.Machine$double.eps)), covtype = "Gaussian")#,
                       #lower = rep(.Machine$double.eps, 5))
homGPmodel <- mleHomGP(X = X, Z = Z, covtype = "Gaussian")
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
test_set_mod$pred <- homGPpred$mean
test_set_mod$se <- sqrt(homGPpred$sd2 + homGPpred$nugs)

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


###for different sites with n_knots = 100
for (site in 1:12) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = site, response = resp, n_knots = 600, n_test_points = 150, package = "hetGP", predictor_type = "normal")
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(site, resp, 600, 150, as.numeric(output[[1]]), as.numeric(output[[2]]), as.numeric(output[[3]]), as.numeric(output[[4]]), "hetGP", "normal")
    #print(model_metrics(site_no = site, response = resp, n_knots = 600, n_test_points = 150, package = "laGP", predictor_type = "normal"))
  }
} 

for (n in c(100, 200, 300, 400, 500, 600)) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = 1, response = resp, n_knots = n, n_test_points = 150, package = "hetGP", predictor_type = "quantile")
    print(output)
    metrics_matrix[nrow(metrics_matrix)+1, ] <- c(1, resp, n, 150, output[[1]], output[[2]], output[[3]], output[[4]], "hetGP", "quantile")
  }
}


for (site in 1:12) {
  for (resp in c(1,2)) {
    output <- model_metrics(site_no = site, response = resp, n_knots = 600, n_test_points = 150, package = "hetGP", predictor_type = "normal", data = shuffled_df)
    metrics_matrix_new[nrow(metrics_matrix_new)+1, ] <- c(site, resp, 600, 150, as.numeric(output[[1]]), as.numeric(output[[2]]), as.numeric(output[[3]]), as.numeric(output[[4]]), "hetGP", "normal")
    #print(model_metrics(site_no = site, response = resp, n_knots = 600, n_test_points = 150, package = "laGP", predictor_type = "normal"))
  }
} 

for (n in c(100, 200, 300, 400, 500, 600)) {
  for (predtype in c("normal", "quantile")) {
    output <- model_metrics(site_no = 1, response = 1, n_knots = n, n_test_points = 150, package = "hetGP", predictor_type = predtype, data = shuffled_df)
    print(output)
    metrics_matrix_new[nrow(metrics_matrix_new)+1, ] <- c(1, resp, n, 150, output[[1]], output[[2]], output[[3]], output[[4]], "hetGP", predtype)
  }
}

save(metrics_matrix_new, file = "metrics_matrix_new.Rdata")
save(metrics_matrix, file = "metrics_matrix.Rdata")

plot1 <- metrics_matrix_new %>%
  filter(package == "hetGP") %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color = knots, shape = predictor_type)) +
  geom_point(size = 3) + labs(x = "time", title = "Effects of increasing knots / changing predictor type", subtitle = "Site No. 1, package used: hetGP")

plot2 <- metrics_matrix_new %>%
  filter(package == "hetGP") %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color = knots, shape = predictor_type)) +
  geom_point(size = 3) + labs(x = "time")

plot11 <- metrics_matrix_new %>%
  filter(response == 1, site == 1, predictor_type == "normal") %>%
  ggplot(mapping = aes(x = time_model + time_pred, y = RMSqE, color = knots, shape = package)) +
  geom_point(size = 3) + labs(x = "time", title = "Differences between packages", subtitle = "site No. = 1, output = 1, predictor_type = normal")

plot12 <- metrics_matrix_new %>%
  filter(response == 1, site == 1, predictor_type == "normal") %>%
  ggplot(mapping = aes(x = time_model + time_pred, y = MStdE, color = knots, shape = package)) +
  geom_point(size = 3) + labs(x = "time")

plot21 <- metrics_matrix%>%
  filter(site == 1, package == "hetGP", response == 1, predictor_type == "normal") %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=RMSqE, color=knots)) +
  geom_point(size = 3) + labs(x = "time", title = "Effects of increasing number of knots using site No.1") #+
  #ylim(0, 1500)

plot22 <- metrics_matrix %>%
  filter(site == 1, package == "hetGP", response == 1, predictor_type == "normal") %>%
  ggplot(mapping = aes(x = time_model+time_pred, y=MStdE, color=knots)) +
  geom_point(size = 3) + labs(x = "time") #+
  #ylim(0,600)


(plot1 <- ggpubr::ggarrange(plot11, plot12, nrow = 2, common.legend = TRUE, legend = "right"))
(plot2 <- ggpubr::ggarrange(plot21, plot22, nrow = 2, common.legend = TRUE, legend = "right"))
(plot <- ggpubr::ggarrange(plot1, plot2, nrow=2, common.legend = TRUE, legend = 'right'))

rmse_normal <- rep(0, 100)
rmse_quantile <- rep(0,100)
mstde_normal <- rep(0,100)
mstde_quantile <- rep(0,100)
time_normal <- rep(0,100)
time_quantile <- rep(0,100)

for (i in 1:100) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_normal <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
                                  package = "hetGP", predictor_type = "normal", data = df)
  metrics_quantile <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150,
                                    package = "hetGP", predictor_type = "quantile", data = df)
  rmse_normal[i] <- metrics_normal[[3]]
  rmse_quantile[i] <- metrics_quantile[[3]]
  mstde_normal[i] <- metrics_normal[[4]]
  mstde_quantile[i] <- metrics_quantile[[4]]
  time_normal[i] <- metrics_normal[[1]] + metrics_normal[[2]]
  time_quantile[i] <- metrics_quantile[[1]] + metrics_quantile[[2]]
}

df_metrics <- data.frame(rmse_normal = rmse_normal, rmse_quantile = rmse_quantile, mstde_normal = mstde_normal,
                         mstde_quantile = mstde_quantile, time_normal = time_normal, time_quantile = time_quantile)


metrics_matrix %>%
  filter(package == "hetGP") %>%
  select(5:8, predictor_type) %>%
  group_by(predictor_type) %>%
  summarise(time = (sum(time_pred) + sum(time_model)) / n() ,
            accuracy = (sum(RMSqE) + sum(MStdE)) / n())


rep <- 50

df_metrics <- data.frame(RMSqE = 1, MStdE = 1, time = 0.1, nan_var = 0.5, nugget = "est")
df_metrics <- df_metrics[-1, ]
df_metrics$nugget <- as.factor(df_metrics$nugget)
levels(df_metrics$nugget) <- c("eps", "0.1", "0.001", "est")

for (i in 1:rep) {
  df <- df[sample(1:nrow(df), size = nrow(df), replace = FALSE), ]
  metrics_eps <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150, package = "hetGP",
                               predictor_type = "normal", nugget_known = TRUE, nugget = sqrt(.Machine$double.eps), data = df)
  metrics_01 <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150, package = "hetGP",
                              predictor_type = "normal", nugget_known = TRUE, nugget = 0.1, data = df)
  metrics_0001 <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150, package = "hetGP",
                                predictor_type = "normal", nugget_known = TRUE, nugget = 0.001, data = df)
  metrics_est <- model_metrics(site_no = 1, response = 1, n_knots = 300, n_test_points = 150, package = "hetGP",
                               predictor_type = "normal", nugget_known = FALSE, data = df)
  
  
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_eps[[3]], metrics_eps[[4]], metrics_eps[[1]]+metrics_eps[[2]],metrics_eps[[5]], "eps")
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_01[[3]], metrics_01[[4]], metrics_01[[1]]+metrics_01[[2]],metrics_01[[5]], "0.1")
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_0001[[3]], metrics_0001[[4]], metrics_0001[[1]]+metrics_0001[[2]],metrics_0001[[5]], "0.001")
  df_metrics[nrow(df_metrics) + 1 , ] <- c(metrics_est[[3]], metrics_est[[4]], metrics_est[[1]]+metrics_est[[2]],metrics_est[[5]], "est")
  
  
}
for (i in 1:4) {
  df_metrics[,i] <- as.numeric(df_metrics[,i])
}

df_metrics %>%
  group_by(nugget) %>%
  summarise(RMSqE = mean(RMSqE),
            MStdE = mean(MStdE, na.rm = TRUE),
            time = mean(time),
            'NaN%' = mean(nan_var)) %>% as.data.frame()
##### EXAMPLES #####


##------------------------------------------------------------
## Example 1: Heteroskedastic GP modeling on the motorcycle data
##------------------------------------------------------------
set.seed(32)

## motorcycle data
library(MASS)
X <- matrix(mcycle$times, ncol = 1)
Z <- mcycle$accel
nvar <- 1
plot(X, Z, ylim = c(-160, 90), ylab = 'acceleration', xlab = "time")

## Model fitting
settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training
model <- mleHetGP(X = X, Z = Z, lower = rep(0.1, nvar), upper = rep(50, nvar),
                  covtype = "Matern5_2", settings = settings)

## A quick view of the fit                  
summary(model)

## Create a prediction grid and obtain predictions
xgrid <- matrix(seq(0, 60, length.out = 301), ncol = 1) 
predictions <- predict(x = xgrid, object =  model)

## Display averaged observations
points(model$X0, model$Z0, pch = 20)

## Display mean predictive surface
lines(xgrid, predictions$mean, col = 'red', lwd = 2)
## Display 95% confidence intervals
lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2)), col = 2, lty = 2)
## Display 95% prediction intervals
lines(xgrid, qnorm(0.05, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
      col = 3, lty = 2)
lines(xgrid, qnorm(0.95, predictions$mean, sqrt(predictions$sd2 + predictions$nugs)), 
      col = 3, lty = 2)

## Comparison with homoskedastic fit
predictions2 <- predict(x = xgrid, object = model$modHom)
lines(xgrid, predictions2$mean, col = 4, lty = 2, lwd = 2)
lines(xgrid, qnorm(0.05, predictions2$mean, sqrt(predictions2$sd2)), col = 4, lty = 3)
lines(xgrid, qnorm(0.95, predictions2$mean, sqrt(predictions2$sd2)), col = 4, lty = 3)



##------------------------------------------------------------
## Example 2: 2D Heteroskedastic GP modeling
##------------------------------------------------------------
set.seed(1)
nvar <- 2

## Branin redefined in [0,1]^2
branin <- function(x){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  x1 <- x[,1] * 15 - 5
  x2 <- x[,2] * 15
  (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
}

## Noise field via standard deviation
noiseFun <- function(x){
  if(is.null(nrow(x)))
    x <- matrix(x, nrow = 1)
  return(1/5*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
}

## data generating function combining mean and noise fields
ftest <- function(x){
  return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
}

## Grid of predictive locations
ngrid <- 51
xgrid <- matrix(seq(0, 1, length.out = ngrid), ncol = 1) 
Xgrid <- as.matrix(expand.grid(xgrid, xgrid))

## Unique (randomly chosen) design locations
n <- 50
Xu <- matrix(runif(n * 2), n)

## Select replication sites randomly
X <- Xu[sample(1:n, 20*n, replace = TRUE),]

## obtain training data response at design locations X
Z <- ftest(X)

## Formating of data for model creation (find replicated observations) 
prdata <- find_reps(X, Z, rescale = FALSE, normalize = FALSE)

## Model fitting
model <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                  lower = rep(0.01, nvar), upper = rep(10, nvar),
                  covtype = "Matern5_2")

model <- mleHetGP(X = X, Z = Z)

## a quick view into the data stored in the "hetGP"-class object
summary(model)                  

## prediction from the fit on the grid     
predictions <- predict(x = Xgrid, object =  model)

## Visualization of the predictive surface
par(mfrow = c(2, 2))
contour(x = xgrid,  y = xgrid, z = matrix(branin(Xgrid), ngrid), 
        main = "Branin function", nlevels = 20)
points(X, col = 'blue', pch = 20)
contour(x = xgrid,  y = xgrid, z = matrix(predictions$mean, ngrid), 
        main = "Predicted mean", nlevels = 20)
points(Xu, col = 'blue', pch = 20)
contour(x = xgrid,  y = xgrid, z = matrix(noiseFun(Xgrid), ngrid), 
        main = "Noise standard deviation function", nlevels = 20)
points(Xu, col = 'blue', pch = 20)
contour(x = xgrid,  y= xgrid, z = matrix(sqrt(predictions$nugs), ngrid), 
        main = "Predicted noise values", nlevels = 20)
points(Xu, col = 'blue', pch = 20)
par(mfrow = c(1, 1))

##------------------------------------------------------------
## Find replicates on the motorcycle data
##------------------------------------------------------------
## motorcycle data
library(MASS)
X <- matrix(mcycle$times, ncol = 1)
Z <- mcycle$accel

data_m <- find_reps(X, Z)

# Initial data
plot(X, Z, ylim = c(-160, 90), ylab = 'acceleration', xlab = "time")
# Display mean values
points(data_m$X0, data_m$Z0, pch = 20)
